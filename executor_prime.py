from argparse import ArgumentParser
from collections import namedtuple
from multiprocessing import Pool
from pickle import dump
from pathlib import Path

Return_tuple    = namedtuple('Return_tuple', ['uniprotkb_ac', 'side_chain_errors', 'backbone_errors', 'pdb2pqr_error_log', 'exit_code'])
Log_tuple       = namedtuple('Log_tuple', ['side_chain_errors', 'backbone_errors', 'pdb2pqr_error_log', 'general_error_log'])


def load_arguments():
    parser = ArgumentParser()
    parser.add_argument('input_PDB_dir',
                        type    = Path,
                        help    = 'Directory with PDB files of structures to be corrected.'
                        )

    parameters = parser.parse_args()
    return parameters

def replica(file: Path):
    try:
        return PrimaryIntegrityMeasuresTaker(input_pdb_file  = file,
                                             from_executor   = True,
                                             log_file        = file.parent/'log.log',
                                             silent          = True,
                                             json_logs_dir   = file.parent/'logs').process_structure()
    except SystemExit as e:
        return Return_tuple(file.stem[3:-12], None, None, None, e.code)
    except:
        return Return_tuple(file.stem[3:-12], None, None, None, 1)

if __name__ == '__main__':
    from prime import PrimaryIntegrityMeasuresTaker

    # load arguments
    args = load_arguments()
    input_dir = args.input_PDB_dir

    # prepare directories
    finished_dir  = input_dir/'finished'
    json_logs_dir = input_dir/'logs'
    finished_dir.mkdir(exist_ok=True)
    json_logs_dir.mkdir(exist_ok=True)

    # prepare data structures
    side_chain_errors   = {}
    backbone_errors     = {}
    pdb2pqr_error_log   = {}
    general_error_log   = []

    # run prime on all files parallelly on 10 CPU cores
    files_n = len([x for x in input_dir.glob('*.pdb')])
    print(f'{files_n}\\0', end='', flush=True)
    with Pool(10) as pool:
        for i, result in enumerate(pool.imap_unordered(replica, input_dir.glob('*.pdb'), chunksize=64), start=1):
            uniprotkb_ac = result.uniprotkb_ac
            if result.side_chain_errors:
                side_chain_errors[uniprotkb_ac] = result.side_chain_errors
            if result.backbone_errors:
                backbone_errors[uniprotkb_ac]   = result.backbone_errors
            if result.pdb2pqr_error_log:
                pdb2pqr_error_log[uniprotkb_ac] = result.pdb2pqr_error_log
            if result.exit_code:
                general_error_log.append((uniprotkb_ac, result.exit_code))

            l = len(str(i-1))
            print(l*'\b' + f'{i}', end='', flush=True)

    # save results
    print('\nINFO: Writing results...', end='')
    with open(args.input_PDB_dir/'results.py', mode='w') as f:
        f.write('side_chain_errors = ')
        f.write(repr(side_chain_errors))
        f.write('\n\nmain_chain_errors = ')
        f.write(repr(backbone_errors))
        f.write('\n\npdb2pqr_error_log = ')
        f.write(repr(pdb2pqr_error_log))
        f.write('\n\ngeneral_error_log = ')
        f.write(repr(general_error_log))
    log_tuple = Log_tuple(side_chain_errors, backbone_errors, pdb2pqr_error_log, general_error_log)
    with open(args.input_PDB_dir/'results.pkl', mode='bw') as f:
        dump(log_tuple, f)
    print('Done.')
