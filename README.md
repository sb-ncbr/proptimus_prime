# PROPTIMUS PRIME

PROPTIMUS PRIME: **P**er-**r**esidue **optim**isation of protein str**u**cture**s**: **Pr**imary **i**ntegrity **me**asures

PROPTIMUS PRIME is a workflow providing repair of amino-acid side chains in protein structures predicted by the [AlphaFold2](https://www.nature.com/articles/s41586-021-03819-2)
algorithm. It can be either run locally as a command line application or integrated as a regular Python library.
It uses the following libraries: [Biopython](https://doi.org/10.1093/bioinformatics/btp163), [RDKit](https://www.rdkit.org/),
 [scikit-learn](https://dl.acm.org/doi/10.5555/1953048.2078195), and the [PDB2PQR](https://doi.org/10.1093/nar/gkh381)
suite.

## Command line use
### Setup
```bash
git clone https://github.com/sb-ncbr/proptimus_prime
cd proptimus_prime/
python -m venv venv_prime
source venv_prime/bin/activate
pip install -r requirements.txt
```

### Execution
#### Single file analysis
```bash
python prime.py <input> [output] [-l log] [-d -s]
```
input: path to the PDB file to process (mandatory)\
output: path to the corrected file (optional)\
log: path to the log file (optional, must be preceded by -l)\
-d: delete auxiliary files (also deletes the log file, if its path has not been specified)\
-s: silent mode

#### Batch run
This runs PRIME over all structures in a folder. This was used for the analysis of Swiss-Prot-proteins predictions in [AlphaFold DB](https://alphafold.com/).
```bash
python executor_prime.py <input_dir> [-c n_cores] 
```
input_dir: path to the directory with PDB files (mandatory)\
n_cores: number of CPU cores to parallelize over

## Python library integration
Import the PrimaryIntegrityMeasuresTaker into your Python script. All options listed in [Command line use](#command-line-use) are
available within the constructor parameters. The logger parameter is used for processing large sets of PDB files using a
custom script. Ignore for casual use.
