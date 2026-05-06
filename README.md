# PROPTIMUS PRIME

PROPTIMUS PRIME: **P**er-**r**esidue **optim**isation of protein str**u**cture**s**: **Pr**imary **i**ntegrity **me**asures

PROPTIMUS PRIME is a workflow providing repair of amino-acid side chains in protein structures predicted by the [AlphaFold2](https://www.nature.com/articles/s41586-021-03819-2)
algorithm. It can be either run locally as a command-line application or integrated as a regular Python library.
It uses the following libraries: [Biopython](https://doi.org/10.1093/bioinformatics/btp163), [RDKit](https://www.rdkit.org/),
 [scikit-learn](https://dl.acm.org/doi/10.5555/1953048.2078195), and the [PDB2PQR](https://doi.org/10.1093/nar/gkh381)
suite.

## Run in a Python virtual environment
### Setup
```bash
git clone https://github.com/sb-ncbr/proptimus_prime
cd proptimus_prime/
python -m venv venv_prime
source venv_prime/bin/activate
pip install -r requirements.txt
```

### Execution
#### Single-file analysis
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

## Running in a Docker container
```bash
# prepare the repository
git clone https://github.com/sb-ncbr/proptimus_prime
cd proptimus_prime/

# build the Docker image
docker build -t local/proptimus-prime .
# or download the image from the registry
docker pull cerit.io/ceitec-biodata-pub/proptimus-prime:latest

# run a single computation in the container
docker run --rm --name proptimus \
  -v ./examples:/opt/proptimus/examples \
  local/proptimus-prime \
  python3 prime.py examples/AF-A4QJE9-F1-model_v6.pdb

# or run the batch run
docker run --rm --name proptimus \
  -v ./examples:/opt/proptimus/examples \
  local/proptimus-prime \
  python3 executor_prime.py examples

# results will be stored in the examples folder
```
## Output
During the [Single-file analysis](), the user is informed about the analysis progress and is provided with the results of the correction.
```
INFO: loading file...
INFO: AF-O15085-F1-model_v4.pdb loaded
INFO: correcting cluster 1: residues 1
INFO: trying with debump...
INFO: success
INFO: correcting cluster 2: residues 895
INFO: success
RESULTS:
Clustered error residues    Correction result
--------------------------  -------------------
1                           success
895                         success
CORRECTION SUCCESSFUL for AF-O15085-F1-model_v4.pdb
INFO: In AF-O15085-F1-model_v4.pdb, backbone errors were found in these residues:
1257
WARNING: Proptimus prime does not provide correction of errors in the backbone.
```
No matter the running mode, information about error detection and correction is written into a log file in the JSON format.

## Python library integration
Import the PrimaryIntegrityMeasuresTaker into your Python script. All options listed in [Execution](#execution) section are
available within the constructor parameters. The logger parameter is used for processing large sets of PDB files using a
custom script. Ignore for casual use.
