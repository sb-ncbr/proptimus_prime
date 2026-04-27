from argparse import ArgumentParser
from collections import defaultdict
from datetime import datetime
from json import dump
from plistlib import InvalidFileException

from networkx import connected_components, Graph
from numpy import array, float32
from numpy.linalg import norm as euclidean_distance
from numpy.typing import NDArray
from os import system
from pathlib import Path
from re import sub
from shutil import rmtree
from tabulate import tabulate

from Bio.PDB import PDBIO, PDBParser, Select
from Bio.PDB.NeighborSearch import NeighborSearch as BioNeighbourSearch
from Bio.PDB.Residue import Residue as BioResidue
from Bio.PDB.Structure import Structure as BioStructure
from rdkit.Chem import MolFromPDBFile
from rdkit.Chem.rdchem import Atom as RDAtom

from data_prime import archetypes, distance, proline_distance, cyclic_residues_atoms

# for batch run use
try:
    from executor_prime import Return_tuple
except ModuleNotFoundError:
    pass
VERSION = '26-04-14'


def print_output(message, silent):
    if not silent:
        print(message)

def save_log(message, log_file):
    with open(log_file, mode='a') as log_file:
        log_file.write(message)


class Atom:
    def __init__(self,
                 rdkit_atom: RDAtom):
        self.rdkit_atom_info = rdkit_atom.GetPDBResidueInfo()
        self.res_id = self.rdkit_atom_info.GetResidueNumber()
        self.res_name = self.rdkit_atom_info.GetResidueName()
        self.name = self.rdkit_atom_info.GetName().strip()
        self.id = self.rdkit_atom_info.GetSerialNumber()
        self.bonded_ats_names: dict[str, Atom] = {}
        self.correct_in_side_chain = True
        self.clashing = False


class Residue:
    def __init__(self,
                 atoms              : dict[int, Atom],
                 missing_bonds      : set[tuple[Atom, str]],
                 side_chain_err_ats : dict[int, Atom]
                 ):

        init_atom                                                   = list(atoms.values())[0]
        self.atoms                                                  = atoms
        self.atoms_by_name          : dict[str, Atom]               = {atom.name: atom for atom in atoms.values()}
        self.clashing                                               = any(atom.clashing for atom in atoms.values())
        self.cyclic_clashing                                        = False
        self.confidence                                             = init_atom.rdkit_atom_info.GetTempFactor()
        self.cycle_centres          : tuple[NDArray[float32], ...]  = (array((0, 0, 0)),)
        self.id                                                     = init_atom.res_id
        self.name                                                   = init_atom.res_name
        self.permeated_cycle                                        = False
        self.side_chain_correct     : bool                          = False if side_chain_err_ats else True
        self.side_chain_err_atoms   : dict[int, Atom]               = side_chain_err_ats if side_chain_err_ats else {}
        self._ox_err_ats            : dict[int, Atom]               = {}

        self.cyclic                                                 = self.name in cyclic_residues_atoms.keys()

        # process missing bond atoms
        self.missing_bonds: dict[tuple[Atom, Atom], NDArray[float32][3]] = {}
        for missing_bond in missing_bonds:
            atom1               = missing_bond[0]
            atom2               = self.atoms_by_name[missing_bond[1]]
            sorted_ats          = sorted((atom1, atom2), key = lambda atom: atom.name)
            sorted_ats_tuple    = (sorted_ats[0], sorted_ats[1])
            self.missing_bonds[sorted_ats_tuple] = array((0, 0, 0))

        # detect oxygen errors
        if side_chain_err_ats:
            # separate backbone oxygen errors (they need to be treated separately)
            for at_id, atom in side_chain_err_ats.items():
                if atom.name in {'O', 'OXT'}:
                    self._ox_err_ats[at_id] = atom
            for at_id in self._ox_err_ats.keys():
                    del self.side_chain_err_atoms[at_id]

            self.min_err_distance = None
            self.calculate_minimal_error_distance()

    def get_kept_ats_ids(self, correction_level:  int) -> set[int]:
        """Returns ids of atoms to be copied into a correction file."""

        # select atoms to be kept
        ats_ids = set()
        # for the case there are only oxygen errors
        if not self.side_chain_err_atoms:
            ats_ids = {atom.id for atom in self.atoms.values() if atom.correct_in_side_chain}
        # for every other case
        else:
            for atom in self.atoms.values():
                # leave out erroneous atoms, including erroneous oxygen atoms
                if atom.correct_in_side_chain:
                    # leave out atoms further in residue than any of the erroneous ones
                    if len(sub(r'[^a-zA-Z]', '', atom.name)) == 2 and atom.name != 'CA':
                        # purvey proline particularities
                        if (self.name == 'PRO'
                                and proline_distance[atom.name[1]] <= self.min_err_distance - correction_level):
                            ats_ids.add(atom.id)
                        elif distance[atom.name[1]] <= self.min_err_distance - correction_level:
                            ats_ids.add(atom.id)
                    else:
                        ats_ids.add(atom.id)

        return ats_ids

    def calculate_minimal_error_distance(self):
        """Calculate "minimal error distance", the distance of the residue's closest side-chain-error atom to the backbone."""
        if self.side_chain_err_atoms:
            # purvey proline particularities
            if self.name == 'PRO':
                self.min_err_distance = min({proline_distance[atom.name[1]]
                                             for atom in self.side_chain_err_atoms.values()})
            else:
                self.min_err_distance = min({distance[atom.name[1]]
                                             for atom in self.side_chain_err_atoms.values()})
        else:
            self.min_err_distance = 1

    def mark_cycle_erroneous(self, cycle_n: int = 1):
        self.side_chain_correct = False
        self.cyclic_clashing    = True
        self.clashing           = True

        for at_name in cyclic_residues_atoms[self.name][cycle_n]:
            atom = self.atoms_by_name[at_name]
            if not self.name == 'PRO' or at_name not in {'N', 'CA'}:
                atom.correct_in_side_chain = False
                self.side_chain_err_atoms[atom.id] = atom

        self.calculate_minimal_error_distance()


class Cluster:
    def __init__(self,
                 err_ress           : set[Residue],
                 inner_ress         : set[Residue],
                 outer_ress         : set[Residue],
                 ):

        self.correct_side_chain_wise    = False
        self.debump                     = False
        self.file                       = ''
        self.inner_ress_ids : set[int]  = {res.id for res in inner_ress}
        self.iterations                 = 0
        self.pdb2pqr_error              = 0
        self.outer_ats_ids  : set[int]  = {atom for res in outer_ress for atom in res.atoms.keys()}

        # determine the side-chain error type (0: other intraresidual, 1: HIS1, 2: HIS2,
        #                                      3: other clash, 4: ring-though clash)
        if any(residue.cyclic_clashing for residue in err_ress):
            self.error_type = 4
        elif any(residue.clashing for residue in err_ress):
            self.error_type = 3
        elif len(err_ress) == 1:
            err_res = list(err_ress)[0]
            if err_res.name == 'HIS':
                pattern = {(atom.name, tuple(at_name for at_name in atom.bonded_ats_names.keys()))
                           for atom in err_res.side_chain_err_atoms.values()}
                if pattern == {('CE1', ()), ('NE2', ())}:
                    self.error_type = 1
                elif pattern == {('CE1', ('ND1',)), ('NE2', ())}:
                    self.error_type = 2
                else:
                    self.error_type = 0
            else:
                self.error_type = 0
        else:
            self.error_type = 3


        sorted_err_ress         : list[Residue]             = sorted(err_ress, key = lambda res: res.id)

        self.err_ress           : dict[int, Residue]        = {res.id: res for res in sorted_err_ress}
        self.confidence_log     : list[tuple[int, float]]   = [(res.id, res.confidence) for res in sorted_err_ress]
        self.err_ress_ids       : list[int]                 = [res.id for res in sorted_err_ress]
        self.err_ress_ids_strs  : list[str]                 = list(map(str, self.err_ress_ids))


class SelectIndexedAtoms(Select):
    """Stores atom ids and uses them within the PDBIO class to select atoms to be copied to the new file
    based on their id."""

    def __init__(self):
        super().__init__()
        self.indices = set()

    def accept_atom(self, atom_id):
        """Overriding the original method (see BioPython documentation)."""
        if atom_id.get_serial_number() in self.indices:
            return 1
        else:
            return 0

    def update_indices(self, indices):
        self.indices = indices


class Protein:
    def __init__(self,
                 path           : Path,
                 corrected_path : Path  = None,
                 correction_dir : Path  = None,
                 cluster_process: bool  = False,
                 final_check    : bool  = False,
                 silent         : bool  = False):
        """The protein is loaded from a PDB file. Error bonds are detected, and the affected atoms and residues are
        noted and clustered."""

        self.filename = path.name
        self.uniprotkb_ac = path.stem[3:-12]
        self._silent = silent

        # load RDKit molecule from PDB file
        rdkit_molecule = MolFromPDBFile(str(path),
                                        sanitize    = False)
        if rdkit_molecule is None:
            raise InvalidFileException(f'ERROR: File at {path} cannot be loaded by RDKit (possibly not a valid PDB file).')

        # make a dictionary of all protein's atoms {atom_id: Atom}, ignoring hydrogen atoms
        atoms : dict[int, Atom] = {rdkit_atom.GetPDBResidueInfo().GetSerialNumber(): Atom(rdkit_atom)
                                   for rdkit_atom in rdkit_molecule.GetAtoms()
                                   if rdkit_atom.GetPDBResidueInfo().GetName().strip()[0] != 'H'}

        # make a set of bonded atoms' NEF names for each atom, meanwhile check for interresidual clashes
        backbone_errs_ress_ids      : set[tuple[int, ...]]  = set()
        side_chain_err_ats          : set[Atom]             = set()
        clashing_ress_pairs_ids     : set[tuple[int, int]]  = set()
        self.ign_proline_ress_ids   : set[int]              = set()
        for bond in rdkit_molecule.GetBonds():

            # ignore hydrogen atoms (added into correction files by propka)
            a1_name = bond.GetBeginAtom().GetPDBResidueInfo().GetName().strip()
            a2_name = bond.GetEndAtom().GetPDBResidueInfo().GetName().strip()
            if a1_name[0] == 'H' or a2_name[0] == 'H':
                continue

            # load atom objects
            atom1 = atoms[bond.GetBeginAtom().GetPDBResidueInfo().GetSerialNumber()]
            atom2 = atoms[bond.GetEndAtom().GetPDBResidueInfo().GetSerialNumber()]

            # purvey disulphide bonds between cysteines' sulfurs
            ats_names = {atom1.name, atom2.name}
            if ats_names == {'SG', 'SG'}:
                continue

            # if the bond is within the same residue or is eupeptidic, add the atom to the set of bonded atoms
            ress_pair = self._make_2_ints_tuple(atom1.res_id, atom2.res_id)
            if atom1.res_id == atom2.res_id:
                if ats_names == {'N', 'C'}:
                    backbone_errs_ress_ids.add((atom1.res_id,))
                    if atom1.res_name == 'PRO':
                        self.ign_proline_ress_ids.add(atom1.res_id)
                else:
                    atom1.bonded_ats_names[atom2.name] = atom2
                    atom2.bonded_ats_names[atom1.name] = atom1
            elif ((atom2.res_id - atom1.res_id == 1 and (atom1.name, atom2.name) == ('C', 'N'))
                    or (atom1.res_id - atom2.res_id == 1 and (atom1.name, atom2.name) == ('N', 'C'))):
                atom1.bonded_ats_names[atom2.name] = atom2
                atom2.bonded_ats_names[atom1.name] = atom1

            # detect interresidual clashes of backbone atoms
            elif ats_names < {'N', 'C', 'CA'}:
                backbone_errs_ress_ids.add(ress_pair)
                if atom1.res_name == 'PRO':
                    self.ign_proline_ress_ids.add(atom1.res_id)

            # otherwise mark atoms as erroneous if they are not backbone atoms
            else:
                both_are_side_chain = True
                for atom in [atom1, atom2]:
                    if atom.name in {'N', 'C', 'CA'}:
                        both_are_side_chain = False
                    else:
                        side_chain_err_ats.add(atom)
                        atom.clashing = True
                        atom.correct_in_side_chain = False

                # mark residues with clashing side chains. the pairs will be used in clustering
                if both_are_side_chain:
                    clashing_ress_pairs_ids.add(ress_pair)

        # check if atoms are bonded to expected atoms
        ress_ids = {atom.res_id for atom in atoms.values()}
        last_res_id = max(ress_ids)
        missing_bonds: set[tuple[int, tuple[Atom, str]]] = set()
        for atom in atoms.values():

            # ignore eventual hydrogen atoms
            at_name = atom.name
            if at_name[0] == 'H':
                continue

            res_name = atom.res_name
            archetype: set[str] = archetypes[res_name][at_name]
            res_id = atom.res_id

            # purvey the terminal carbon and oxygen
            if res_id == last_res_id:
                if at_name == 'C':
                    archetype = archetype - {'N'}
                    archetype.add('OXT')
                if at_name == 'OXT':
                    archetype = archetypes[res_name]['O']

            # purvey the initial hydrogen and the ends of local chains in correction files
            elif (at_name == 'N' and res_id-1 not in ress_ids
                  or at_name == 'C' and res_id+1 not in ress_ids):
                archetype = archetype - {'N', 'C'}

            # select atoms to be marked as erroneous
            bonded_ats_names = atom.bonded_ats_names.keys()
            correct_in_side_chain = True
            if bonded_ats_names != archetype:
                # if backbone oxygen, mark as erroneous
                if at_name in {'O', 'OXT'}:
                    correct_in_side_chain = False

                else:
                    # note backbone errors
                    if at_name == 'CA' and not {'N', 'C'} <= bonded_ats_names:
                        backbone_errs_ress_ids.add((res_id,))
                        if res_name == 'PRO' and 'N' not in bonded_ats_names:
                            self.ign_proline_ress_ids.add(res_id)

                    elif at_name == 'N' and {'C'} & archetype != {'C'} & bonded_ats_names:
                        backbone_errs_ress_ids.add((res_id - 1, res_id))

                    # find side-chain errors, thus ignore backbone atoms
                    elif at_name not in {'N', 'C', 'CA'}:

                        # if there are any extra atoms bonded, mark erroneous - intraresidual clashes detection
                        if bonded_ats_names - archetype:
                            correct_in_side_chain = False

                        # purvey proline particularities
                        if res_name == 'PRO':
                            if at_name == 'CD' and 'N' not in bonded_ats_names:
                                correct_in_side_chain = False
                            elif at_name == 'CG':
                                correct_in_side_chain = False

                        # ignore the atoms that are the "very last correct" in the residue
                        elif not all(distance[missing_at_nef_name[1]] > distance[at_name[1]]
                                     for missing_at_nef_name in archetype - bonded_ats_names
                                     if len(sub(r'[^a-zA-Z]', '', missing_at_nef_name)) == 2): # remove an eventual digit at the end of the nef name (e.g. CE1 and NE2 in HIS)
                            correct_in_side_chain = False

                if not correct_in_side_chain:
                    atom.correct_in_side_chain = False
                    side_chain_err_ats.add(atom)

                # note missing bonds
                for unbonded_at_name in archetype - bonded_ats_names:
                    if (at_name, unbonded_at_name) != ('N', 'C'): # note missing eupeptide bonds only once
                        if (at_name, unbonded_at_name) == ('C', 'N'):
                            missing_bonds.add((res_id + 1, (atom, unbonded_at_name)))
                        else:
                            missing_bonds.add((res_id, (atom, unbonded_at_name)))

        # summarize heretofore results
        self.side_chains_correct    = False if side_chain_err_ats else True
        self.backbone_correct       = False if backbone_errs_ress_ids else True
        ress_ats = self._make_ress_ats_dicts_dict(set(atoms.values()))
        side_chain_err_ress_ats = self._make_ress_ats_dicts_dict(side_chain_err_ats)
        self.residues: dict[int, Residue]   = {res_id: Residue(atoms              = res_ats_dict,
                                                               missing_bonds      = {bond[1] for bond in missing_bonds
                                                                                     if bond[0] == res_id},
                                                               side_chain_err_ats = side_chain_err_ress_ats[res_id] if res_id in side_chain_err_ress_ats.keys()
                                                                                    else dict())
                                               for res_id, res_ats_dict in ress_ats.items()}
        self.clusters: list[Cluster]        = []

        # prepare for correction if necessary
        if not cluster_process: # supposing pdb2pqr would not create ring-through clashes, we do not need to search for them in cluster processes
            self.side_chain_err_ress_ids            : set[int] = set()
            if not self.side_chains_correct or not self.backbone_correct:
                self._bio_structure     : BioStructure                              = PDBParser(QUIET=True).get_structure('protein', path)
                cycle_distance                                                      = 1
                self.log                                                            = ''
                self.pdb2pqr_errors_log : list[tuple[int, list[tuple[int, float]]]] = []
                self._surroundings_distance                                         = 15

                # find neighbouring residues to the erroneous ones and find errors related to aromatic rings
                self.side_chain_err_ress_ids = set(side_chain_err_ress_ats.keys())
                del side_chain_err_ress_ats
                self._bio_residues  : dict[int, BioResidue]             = {residue.id[1]: residue
                                                                           for residue in self._bio_structure.get_residues()}
                self._kdtree        : BioNeighbourSearch                = BioNeighbourSearch(list(self._bio_structure.get_atoms()))
                neighbouring_ress   : dict[int, tuple[set[Residue],
                                                      set[Residue]]]    = {}
                for residue in self.residues.values():

                    # find neighbours for every residue with a side chain error
                    if not residue.side_chain_correct:
                        neighbouring_ress[residue.id] = (self._find_neighbours_ids(subject = residue,
                                                                                   find_close_neighbours = True),
                                                         self._find_neighbours_ids(subject = residue))

                    # ignore residues not missing any bond
                    if not residue.missing_bonds:
                        continue

                    # determine centres of overextended bonds
                    for at_pair in residue.missing_bonds.keys():
                        bio_bond_at1 = self._bio_structure[0]['A'][(' ', at_pair[0].res_id, ' ')][at_pair[0].name]
                        bio_bond_at2 = self._bio_structure[0]['A'][(' ', at_pair[1].res_id, ' ')][at_pair[1].name]
                        x = (bio_bond_at1.coord[0] + bio_bond_at2.coord[0])/2
                        y = (bio_bond_at1.coord[1] + bio_bond_at2.coord[1])/2
                        z = (bio_bond_at1.coord[2] + bio_bond_at2.coord[2])/2
                        residue.missing_bonds[at_pair] = array((x, y, z))

                    # find the closest neighbouring residues among which eventual cyclic clash partners will be searched
                    if residue.id not in neighbouring_ress.keys():
                        close_neigh_ress = self._find_neighbours_ids(subject            = residue,
                                                                     find_close_neighbours= True)
                    else:
                        close_neigh_ress = neighbouring_ress[residue.id][0]
                    close_neigh_cycl_ress = {res for res in close_neigh_ress if res.cyclic}
                    if not close_neigh_cycl_ress:
                        continue

                    # determine geometric centres of cycles in the closest neighbouring cyclic residues
                    for neigh_cycl_res in close_neigh_cycl_ress:
                        centres: list[NDArray[float32]] = []
                        for cycle in cyclic_residues_atoms[neigh_cycl_res.name]:
                            bio_cycle_atoms = [self._bio_structure[0]['A'][(' ', neigh_cycl_res.id, ' ')][at_name]
                                               for at_name in cycle]
                            x = sum(atom.coord[0] for atom in bio_cycle_atoms)/len(bio_cycle_atoms)
                            y = sum(atom.coord[1] for atom in bio_cycle_atoms)/len(bio_cycle_atoms)
                            z = sum(atom.coord[2] for atom in bio_cycle_atoms)/len(bio_cycle_atoms)
                            centres.append(array((x, y, z)))
                        neigh_cycl_res.cycle_centres = tuple(centre for centre in centres)

                    # detect the cyclic clashes based on coordinates
                    for neigh_cycl_res in close_neigh_cycl_ress:
                        for at_pair, centre in residue.missing_bonds.items():
                            for cycle_id, cycle_centre in enumerate(neigh_cycl_res.cycle_centres):

                                # decide whether the overextended bond goes through the cycle
                                if euclidean_distance(centre - cycle_centre) < cycle_distance:
                                    residue.cyclic_clashing = True
                                    ress_pairs = [self._make_2_ints_tuple(residue.id, neigh_cycl_res.id)]

                                    # purvey for ring deformation in tyrosine due to CZ-OH bond going through a ring
                                    if 'OH' in (atom.name for atom in at_pair): # the OH atom is only found in tyrosine
                                        residue.mark_cycle_erroneous()

                                    # purvey for cases of two rings clashing
                                    elif residue.cyclic:
                                        for c_id, cycle_ats in enumerate(cyclic_residues_atoms[residue.name]):
                                            if {at.name for at in at_pair} < set(cycle_ats):
                                                residue.mark_cycle_erroneous(c_id)
                                                residue.permeated_cycle = True

                                    # assure proline-through clashes be processed only on the level of the other residue
                                    if neigh_cycl_res.name == 'PRO':
                                        self.ign_proline_ress_ids.add(neigh_cycl_res.id)
                                        pair_names = {a.name for a in at_pair}

                                        # purvey the case of purely backbone errors
                                        if (residue.name == 'PRO' and 'O' not in [name[0] for name in pair_names]
                                                or pair_names < {'N', 'CA', 'C'}):
                                            ress_ids = ress_pairs[0]

                                            # purvey the case of 2 residues separated by proline
                                            if pair_names == {'C', 'N'}:
                                                pair_ress_ids = [at.res_id for at in at_pair]
                                                ress_ids = self._make_ints_tuple([*pair_ress_ids] + [neigh_cycl_res.id])
                                                backbone_errs_ress_ids.add(ress_ids)
                                                backbone_errs_ress_ids.discard((pair_ress_ids[0], pair_ress_ids[1]))
                                                backbone_errs_ress_ids.discard((pair_ress_ids[1], pair_ress_ids[0]))
                                                ress_pairs.extend({self._make_2_ints_tuple(pair_res_id, neigh_cycl_res.id)
                                                                   for pair_res_id in pair_ress_ids})
                                            else:
                                                backbone_errs_ress_ids.add(ress_ids)
                                            for res_id in ress_ids:
                                                backbone_errs_ress_ids.discard((res_id,))
                                        else:
                                            backbone_errs_ress_ids.add((neigh_cycl_res.id,))
                                        break

                                    # apply cyclic residue as an erroneous residue
                                    neigh_cycl_res.mark_cycle_erroneous(cycle_id)
                                    neigh_cycl_res.permeated_cycle = True
                                    if neigh_cycl_res.id not in self.side_chain_err_ress_ids:
                                        self.side_chain_err_ress_ids.add(neigh_cycl_res.id)
                                        neighbouring_ress[neigh_cycl_res.id] = (
                                            self._find_neighbours_ids(subject               = neigh_cycl_res,
                                                                      find_close_neighbours = True),
                                            self._find_neighbours_ids(subject               = neigh_cycl_res))

                                    # add to clashing residues pairs
                                    clashing_ress_pairs_ids.add(*ress_pairs)

                                    break # even if the bond would go through both rings of W, nothing would change

                # update protein status and neighbours data
                self.side_chain_err_ress_ids = {res_id for res_id in self.side_chain_err_ress_ids if res_id not in self.ign_proline_ress_ids}
                self.side_chains_correct = False if self.side_chain_err_ress_ids else True
                side_chain_or_ign_pro_ress = {self.residues[res_id] for res_id in self.side_chain_err_ress_ids | self.ign_proline_ress_ids}
                for res_id, neigh_ress in neighbouring_ress.items():
                    neighbouring_ress[res_id] = (neigh_ress[0] - side_chain_or_ign_pro_ress,
                                                 neigh_ress[1])

                # make clusters
                if len(self.side_chain_err_ress_ids) == 1:
                    err_res = self.residues[list(self.side_chain_err_ress_ids)[0]]
                    self.clusters = [Cluster({err_res},
                                             neighbouring_ress[err_res.id][0] | {err_res},
                                             neighbouring_ress[err_res.id][1] - {err_res})]
                else:
                    edges = []
                    for res_id_1, res_id_2 in clashing_ress_pairs_ids:
                        if {res_id_1, res_id_2} <= self.side_chain_err_ress_ids:
                            edges.append((res_id_1, res_id_2))
                    graph = Graph()
                    graph.add_edges_from(edges)
                    graph.add_nodes_from(self.side_chain_err_ress_ids)
                    for cluster in list(connected_components(graph)):
                        cluster_err_ress                    = {self.residues[res_id]
                                                               for res_id in cluster}
                        close_neigh_ress_sets               = [neighbouring_ress[res_id][0] for res_id in cluster]
                        neigh_ress_sets                     = [neighbouring_ress[res_id][1] for res_id in cluster]
                        cluster_inner_ress  : set[Residue]  = set().union(*close_neigh_ress_sets) | cluster_err_ress
                        cluster_outer_ress  : set[Residue]  = set().union(*neigh_ress_sets) - cluster_err_ress
                        self.clusters.append(Cluster(cluster_err_ress, cluster_inner_ress, cluster_outer_ress))

            # summarize heretofore results
            self.backbone_correct                                       = False if backbone_errs_ress_ids else True
            self.backbone_errors        : list[tuple[int, ...]]         = sorted(backbone_errs_ress_ids, key = lambda res: res[0])
            self.backbone_err_ress_ids  : list[int]                     = sorted({res_id
                                                                                  for err in backbone_errs_ress_ids
                                                                                  for res_id in err})
            if not final_check:
                if not correction_dir or not corrected_path:
                    raise NotImplementedError('The main protein object must have its correction_dir and corrected_path specified.')
                self._correction_dir    : Path                          = correction_dir if correction_dir else Path()
                self._final_file_path   : Path                          = corrected_path

    @staticmethod
    def _make_ress_ats_dicts_dict(atoms: set[Atom]) -> dict[int, dict[int, Atom]]:
        ress_ats_dicts_dict = defaultdict(dict[int, Atom])
        for atom in atoms:
            ress_ats_dicts_dict[atom.res_id][atom.id] = atom
        return dict(ress_ats_dicts_dict)

    @staticmethod
    def _make_2_ints_tuple(member1: int, member2: int) -> tuple[int, int]:
        sorted_ints = sorted([member1, member2])
        return sorted_ints[0], sorted_ints[1]

    @staticmethod
    def _make_ints_tuple(members: list[int]) -> tuple[int, ...]:
        sorted_ints = sorted(members)
        return tuple(sorted_ints)

    def _find_neighbours_ids(self,
                             subject                : Residue,
                             find_close_neighbours  : bool = False) -> set[Residue]:
        neighbours : set[Residue]  = set()
        distance_divisor = 2 if find_close_neighbours else 1
        for residue in set(self._kdtree.search(center = self._bio_residues[subject.id].center_of_mass(geometric=True),
                                               radius = self._surroundings_distance / distance_divisor,
                                               level  = 'R')):
            res = self.residues[residue.id[1]]
            if res != subject:
                neighbours.add(res)

        return neighbours

    def execute_correction(self):
        """Execute correction over individual clusters"""

        # try iterations to correct the clusters
        io = PDBIO()
        io.set_structure(self._bio_structure)
        selector = SelectIndexedAtoms()
        for cluster_id, cluster in enumerate(self.clusters, start=1):
            print_output(f'INFO: correcting cluster {cluster_id}: residues {'+'.join(cluster.err_ress_ids_strs)}', self._silent)

            # ensure a directory exists for the current correction level
            cluster_corr_dir = self._correction_dir/f'cluster_{cluster_id}'
            cluster_corr_dir.mkdir(exist_ok=True)

            # try iterations of correction
            max_error_distance = max({res.min_err_distance if res.min_err_distance else 1
                                      for res in cluster.err_ress.values()})
            correction_level = 1
            correction_attempt = self
            output_file = Path()
            pdb2pqr_problem = False
            pdb2pqr_nodebump = '--nodebump'
            debump_flag = ''
            while not all(correction_attempt.residues[res_id].side_chain_correct
                          for res_id in cluster.inner_ress_ids):
                if correction_level > max_error_distance or pdb2pqr_problem:
                    if not pdb2pqr_problem and pdb2pqr_nodebump:
                        print_output('INFO: trying with debump...', self._silent)
                        correction_level = 1
                        pdb2pqr_nodebump = ''
                        debump_flag = '_d'
                        cluster.debump = True
                    else:
                        print_output(f'INFO: failure', self._silent)
                        cluster.iterations = correction_level - 1
                        break
                cutout_file = cluster_corr_dir/f'level{correction_level}.pdb'
                output_file = cluster_corr_dir/f'level{correction_level}{debump_flag}_out.pdb'
                pdb2pqr_log = cluster_corr_dir/f'level{correction_level}{debump_flag}_pdb2pqr.log'

                if pdb2pqr_nodebump:
                    # cut out the cluster's atoms into a correction file
                    kept_ats_ids = set()
                    for res_id in cluster.err_ress.keys():
                        kept_ats_ids.update(self.residues[res_id].get_kept_ats_ids(correction_level))
                    surr_ats_ids = cluster.outer_ats_ids
                    selector.update_indices(kept_ats_ids | surr_ats_ids)
                    io.save(str(cutout_file), selector)

                # correct by propka
                system(f'pdb2pqr30 {pdb2pqr_nodebump} --noopt --pdb-output {output_file} '
                       f'{cutout_file} {cluster_corr_dir}/delete.pqr '
                       f'2>{pdb2pqr_log};'
                       f'rm {cluster_corr_dir}/delete.*;')

                # check if the outcut cluster was not too small for pdb2pqr
                with open(pdb2pqr_log, mode='r') as log_file:
                    for line in log_file:
                        if line[0:34] == 'ERROR:This PDB file is missing too':
                            self.pdb2pqr_errors_log.append((cluster_id, cluster.confidence_log))
                            cut = line.find('T')
                            number = line.find('0.')
                            num_end = line.find(')')
                            cluster.pdb2pqr_error = float(line[number:num_end])
                            print_output(f'ERROR: The level {correction_level} cutout {line[15:cut]} {line[cut:-1]}', self._silent)
                            self.log += f'ERROR: The level {correction_level} cutout {line[15:108]} {line[108:]}'

                            pdb2pqr_problem = True
                            cluster.iterations = correction_level
                            break

                # load the correction attempt output file to check successfulness
                if not pdb2pqr_problem:
                    correction_attempt = Protein(path           = output_file,
                                                 cluster_process= True)
                correction_level += 1

            else:
                print_output('INFO: success', self._silent)
                for res_id in cluster.err_ress.keys():
                    self.residues[res_id].side_chain_correct = True

                cluster.correct_side_chain_wise = True
                cluster.iterations = correction_level - 1
                cluster.file = f'{cluster_corr_dir}{debump_flag}.pdb'
                output_file.replace(cluster.file)

        # check if all side chain errors were corrected
        if all(cluster.correct_side_chain_wise for cluster in self.clusters):
            self.side_chains_correct = True

        # write the corrected cluster into the output file
        for i, cluster in enumerate((cluster for cluster in self.clusters if cluster.correct_side_chain_wise), start=1):
            cluster_bio : BioStructure = PDBParser(QUIET=True).get_structure(id    = f'cluster{i}',
                                                                             file  = cluster.file)
            for atom in cluster_bio.get_atoms():
                res_id = atom.get_parent().id[1]
                # write only the inner parts of the cluster
                if res_id in cluster.inner_ress_ids:
                    # remove extra atoms created by pdb2pqr (OXT atom)
                    if atom.name in self._bio_structure[0]['A'][res_id]:
                        self._bio_structure[0]['A'][res_id][atom.name].coord = atom.coord
        io.save(str(self._final_file_path))


class PrimaryIntegrityMeasuresTaker:
    def __init__(self,
                 input_pdb_file         : Path,
                 output_pdb_file        : Path              = None,
                 from_executor                              = False,
                 log_file               : Path              = None,
                 delete_auxiliary_files                     = False,
                 silent                                     = True,
                 json_logs_dir          : Path              = None):

        # control usability of files
        if not input_pdb_file.is_file():
            raise FileNotFoundError(f'ERROR: File at {input_pdb_file} does not exist!')
        if input_pdb_file == output_pdb_file:
            raise ValueError(f'ERROR: Input and output files cannot be the same.')

        self._input_pdb_file = input_pdb_file
        prime_dir = input_pdb_file.parent/'correction_prime'
        self._correction_dir = prime_dir/(input_pdb_file.stem[3:-12]+'_correction')
        if output_pdb_file:
            self._output_PDB_file = output_pdb_file
        else:
            self._output_PDB_file = prime_dir/(input_pdb_file.stem + '_corrected.pdb')
        self._output_PDB_file.parent.mkdir(exist_ok=True, parents=True)

        if json_logs_dir:
            self.json_logs_dir = json_logs_dir
        else:
            self.json_logs_dir = prime_dir
        self.json_logs_dir.mkdir(parents = True, exist_ok = True)

        self._from_executor = from_executor
        self._log = ''
        self._log_file = log_file
        self._delete_auxiliary_files = delete_auxiliary_files
        self._silent = silent

    def process_structure(self) -> Return_tuple:
        # load protein into a Protein object, check for errors
        print_output('INFO: loading file...', self._silent)
        protein = Protein(path              = self._input_pdb_file,
                          corrected_path    = self._output_PDB_file,
                          correction_dir    = self._correction_dir,
                          silent            = self._silent)
        print_output(f'INFO: {protein.filename} loaded', self._silent)

        side_chain_errors   = ()
        pdb2pqr_error_log   = []
        backbone_errors     = []
        erroneous_correction = False
        if protein.side_chains_correct and protein.backbone_correct:
            print_output(f'OK: No error found in {protein.filename}.', self._silent)
        else:
            self._log += f'{protein.uniprotkb_ac}:\n'
            if not protein.side_chains_correct:
                self._correction_dir.mkdir(exist_ok=True, parents=True)
                protein.execute_correction()

                # run a check of the final file
                persistent_side_chain_err_ress_ids = (protein.side_chain_err_ress_ids
                                                      - {res_id
                                                         for cluster in protein.clusters if cluster.correct_side_chain_wise
                                                         for res_id in cluster.err_ress.keys()})
                final_protein = Protein(path        = self._output_PDB_file,
                                        final_check = True)
                side_chain_errors_different = final_protein.side_chain_err_ress_ids != persistent_side_chain_err_ress_ids
                backbone_errors_different = set(final_protein.backbone_errors) != set(protein.backbone_errors)
                if side_chain_errors_different or backbone_errors_different:

                    # purvey ignored prolines classified differently in final_protein because the other residue's
                    # clashing side chain was corrected
                    erroneous_correction = False
                    new_backbone_errors = {new_err for new_err in final_protein.backbone_errors
                                           if not any(set(new_err) <= set(old_err) for old_err in protein.backbone_errors)}
                    if new_backbone_errors:
                        erroneous_correction = True
                    else:
                        new_side_chain_err_ress_ids = final_protein.side_chain_err_ress_ids - persistent_side_chain_err_ress_ids
                        disapp_backbone_err_pro_ids = protein.ign_proline_ress_ids - set(final_protein.backbone_err_ress_ids)
                        supposed_benign_new_side_chain_err_ress_ids = new_side_chain_err_ress_ids & disapp_backbone_err_pro_ids
                        supposed_benign_but_a_new_clash = any(final_protein.residues[res_id].clashing
                                                              for res_id in supposed_benign_new_side_chain_err_ress_ids)
                        surely_malign_new_side_chain_errs = (new_side_chain_err_ress_ids - supposed_benign_new_side_chain_err_ress_ids) != set()
                        malign_new_side_chain_errs = surely_malign_new_side_chain_errs or supposed_benign_but_a_new_clash
                        if malign_new_side_chain_errs:
                            erroneous_correction = True

                    if erroneous_correction:
                        print_output('ERROR: pdb2pqr failed to correct the protein.', self._silent)
                        for cluster in protein.clusters:
                            cluster.correct_side_chain_wise = False
                        protein.side_chains_correct = False
                        erroneous_correction = True

                # log the results
                self._log += protein.log
                listed_results = [[', '.join(cluster.err_ress_ids_strs), 'success' if cluster.correct_side_chain_wise else 'failure']
                                  for cluster in protein.clusters]
                table = tabulate(tabular_data   = listed_results,
                                 headers        = ['Clustered error residues', 'Correction result'],
                                 colalign       = ['left', 'left'])
                if not erroneous_correction:
                    print_output(f'RESULTS:\n{table}', self._silent)
                self._log += f'{table}\n'
                if protein.side_chains_correct:
                    print_output(f'CORRECTION SUCCESSFUL for {protein.filename}', self._silent)
                elif all(not cluster.correct_side_chain_wise for cluster in protein.clusters):
                    print_output(f'CORRECTION FAILED for {protein.filename}', self._silent)
                else:
                    print_output(f'CORRECTION partially successful for {protein.filename}', self._silent)

                if self._from_executor:
                    side_chain_errors = ([(cluster.confidence_log,
                                           cluster.error_type,
                                           cluster.correct_side_chain_wise,
                                           cluster.iterations,
                                           cluster.debump,
                                           cluster.pdb2pqr_error,
                                           )
                                          for cluster in protein.clusters],
                                         protein.side_chains_correct)

                    if protein.pdb2pqr_errors_log:
                        pdb2pqr_error_log = protein.pdb2pqr_errors_log

            if not protein.backbone_correct:
                print_output(f'INFO: In {protein.filename}, backbone errors were found in these residues:', self._silent)
                chain_errors_string = ', '.join(list(map(str, protein.backbone_err_ress_ids)))
                print_output(chain_errors_string, self._silent)
                self._log += (f'--------------------------\n'
                              f'Chain errors\n'
                              f'--------------------------\n'
                              f'{chain_errors_string}\n')
                print_output('WARNING: Proptimus prime does not provide correction of errors in the backbone.', self._silent)

                if self._from_executor:
                    backbone_errors = [[(res_id, protein.residues[res_id].confidence) for res_id in error]
                                       for error in protein.backbone_errors]

            if erroneous_correction:
                self._log += 'ERROR: PDB2PQR failed to correct the protein.\n'

            # write the log or delete auxiliary files (if demanded)
            self._log += '\n'
            if self._delete_auxiliary_files:
                if self._correction_dir.is_dir():
                    rmtree(self._correction_dir)
            else:
                if self._log_file:
                    self._log_file.parent.mkdir(exist_ok = True, parents = True)
                else:
                    self._log_file = self._correction_dir/'log.txt'
                    self._log_file.parent.mkdir(exist_ok = True, parents = True)
                with open(self._log_file, mode='a') as log_file:
                    log_file.write(self._log)

        # log into a .json file
        if not self._from_executor:
            json_log = {'file': str(self._input_pdb_file),
                        'date_and_time': str(datetime.now()),
                        'version': VERSION,
                        'side_chain_errors':
                            {'detected_n': len(protein.clusters),
                             'repaired_n': len([cluster for cluster in protein.clusters
                                                if cluster.correct_side_chain_wise]),
                             'list': [{'affected_residues': cluster.err_ress_ids,
                                       'repaired': cluster.correct_side_chain_wise}
                                      for cluster in protein.clusters]},
                        'backbone_errors':
                            {'detected_n': len(protein.backbone_errors),
                             'list': [{'affected_residues': backbone_error}
                                      for backbone_error in protein.backbone_errors]}}
        else:
            json_log = {'file': str(self._input_pdb_file),
                        'date_and_time': str(datetime.now()),
                        'version': VERSION,
                        'erroneous correction': erroneous_correction,
                        'side_chain_errors':
                            {'detected_n': len(protein.clusters),
                             'repaired_n': len([cluster for cluster in protein.clusters
                                                if cluster.correct_side_chain_wise]),
                             'list': [{'affected_residues': [{'residue_id': tup[0],
                                                              'confidence': tup[1]}
                                                             for tup in cluster.confidence_log],
                                       'repaired': cluster.correct_side_chain_wise,
                                       'error_type': cluster.error_type,
                                       'iterations_n': cluster.iterations,
                                       'debump': cluster.debump,
                                       'pdb2pqr_error': cluster.pdb2pqr_error}
                                      for cluster in protein.clusters]},
                        'backbone_errors':
                            {'detected_n': len(protein.backbone_errors),
                             'list': [{'affected_residues': backbone_error}
                                      for backbone_error in protein.backbone_errors]}}

        with open(self.json_logs_dir/(protein.uniprotkb_ac + '_log.json'), 'w', encoding='utf8') as f:
            dump(json_log, f)

        if self._from_executor:
            return Return_tuple(protein.uniprotkb_ac, side_chain_errors, backbone_errors, pdb2pqr_error_log, 5 if erroneous_correction else 0)

        if erroneous_correction:
            raise RuntimeError
        else:
            return Return_tuple(0, 0, 0, 0, 0)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('input_PDB_file',
                        type    = Path,
                        help    = 'PDB file with structure to be corrected.'
                        )
    parser.add_argument('output_PDB_file',
                        type    = Path,
                        nargs   = '?',
                        default = None,
                        help    = 'Path for the corrected structure file. (optional)'
                        )
    parser.add_argument('-l', '--log_file',
                        type    = Path,
                        nargs   = '?',
                        default = None,
                        help    = 'Path for the logging file.')
    parser.add_argument('-d', '--delete_auxiliary_files',
                        action  = 'store_true',
                        help    = 'Delete all auxiliary files. (recommended)'
                        )
    parser.add_argument('-s', '--silent',
                        action  = 'store_true',
                        help    = 'Reduce output. (use only if you know what you are doing)'
                        )

    args = parser.parse_args()
    PrimaryIntegrityMeasuresTaker(args.input_PDB_file,
                                  args.output_PDB_file,
                                  False,
                                  args.log_file,
                                  args.delete_auxiliary_files,
                                  args.silent).process_structure()
