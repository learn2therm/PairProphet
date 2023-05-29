"""
Main script wrapper.
This is pseudocode for our parent script.

Directory: scipts

Input: 

Output: 

Note: 
We have two wrapping functions, but 
these can be developed as individual scripts.

Runtime: 
"""
# system dependencies
import sys
import logging
import os

# library dependencies
import pandas as pd
from joblib import Parallel, delayed

# local dependencies

##machine learning
from PairPro.evaluate_model import evaluate_model
from PairPro.train_val_wrapper import train_val_wrapper
from PairPro.train_val_input_cleaning import columns_to_keep

##build DB
from PairPro.preprocessing import connect_db, build_fafsa
from PairPro.user_blast import make_blast_df

##hmmer
import PairPro.hmmer
import PairPro.utils


##structure
# from PairPro.structure import download_structures, run_fatcat


## HMMER Paths
HMM_PATH = './data/pfam/Pfam-A.hmm'  # ./Pfam-A.hmm
PRESS_PATH = './data/pfam'
HMMER_OUTPUT_DIR = './data/protein_pairs/'
WORKER_WAKE_UP_TIME = 25 # this is to ensure that if a worker that is about to be shut down due to previous task completetion doesn't actually start running

## get environmental variables
if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
    LOGLEVEL = getattr(logging, LOGLEVEL)
else:
    LOGLEVEL = logging.INFO
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'

#Ryans Component
"""
Input: learn2thermDB (queuing format)
Output: SQL db after filtering, and then each component can retrieve the data as needed.
Params: Base environment + DuckDB
"""
con, _ = connect_db(dbpath)
build_fafsa(con)

#Humood/Amin Component
"""
Input: SQL db from Ryan, then retrieve only 1000 pairs
Output: SQL db with appended Boolean from running HMMR remotely
Params: Base environment + PyHMMER + Joblib
"""
df = con.execute("""SELECT m_protein_seq, t_protein_seq, prot_pair_index, meso_pid, thermo_pid, meso_pdb, thermo_pdb""").df()
# num_pids = con.execute("""SELECT COUNT(pid) FROM pairpro_proteins""").df()

# # get the number of sequences
# nseqs = len(num_pids)
# if nseqs < 1000:
#     logger.info(f'Number of sequences is less than 1000, so all sequences will be processed via API')
#     # run hmmscan
#     PairPro.hmmer.hmmerscanner()
#     PairPro.hmmer.run_hmmerscanner()
# else:
# logger.info(f'Number of sequences is greater than 1000, so sequences will be processed in chunks locally')
proteins_in_pair_pids = con.execute("SELECT pid FROM fafsa_proteins").df()
logger.debug(f"Total number of protein in pairs: {len(proteins_in_pair_pids)} in pipeline")

# chunking the PID so the worker function queries
protein_pair_pid_chunks = [proteins_in_pair_pids[i:i + chunk_size]
                            for i in range(0, len(proteins_in_pair_pids), chunk_size)]

# run hmmscan
logger.info('Running pyhmmer in parallel on all chunks')

Parallel(
    n_jobs=njobs)(
    delayed(PairPro.hmmer.local_hmmer_wrapper)(
        chunk_index,
        dbpath,
        protein_pair_pid_chunks,
        None) for chunk_index,
    protein_pair_pid_chunks in enumerate(protein_pair_pid_chunks))

logger.debug(f"number of protein chunks: {len(protein_pair_pid_chunks)}")

#Structure Component
"""
Input: SQL db from Ryan, retrieve only PDB IDs and Uniprot IDs (sequences are not required for this component)
Output: SQL db appended with Boolean on whether the structures are similar
Params: Base environment + FATCAT
"""
# download_structures(df)
# run_fatcat(df)

#Logan Component
"""
Input: SQL db from Humood or Amin (with HMMR Boolean + prot_pair_ind + Jaccard score appended)
Output: csv files with 6 columns (seq1, seq2, protein_match (Humood + Amin), protein_match(ML), Jaccard Score)
Params: Base environment + iFeatureOmega dependencies 
"""
#need to save model to use in next script
model_output = train_val_wrapper(df, feature_list=[])
con.commit()
con.close()

trained_model = model_dev(dbpath)[2]


if __name__ == "__main__":
    # Initialize logger
    logger = PairPro.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')
    logger.info(f"Running {__file__}")

    # create pfam HMM directory (this was before HMM download script)
    try:
        os.makedirs('./data/pfam', exist_ok=True)
    except OSError as e:
        logger.error(f'Error creating directory: {e}')

    # press the HMM db
    PairPro.hmmer.hmmpress_hmms(HMM_PATH, PRESS_PATH)

    logger.info(f'Pressed HMM DB: {PRESS_PATH}')

    dbpath = "TO DO"

    #connect to database
    con, _ = connect_db(dbpath)
    build_fafsa(con)
    logger.info('Connected to database.')

    # if local hmmer
    # Set up parallel processing and parsing
    chunk_size = int(sys.argv[1]) # Number of sequences to process in each chunk
    njobs = int(sys.argv[2])  # Number of parallel processes to use

    logger.info('Parallel processing parameters obtained')

    # prepare output file
    try:
        os.makedirs(HMMER_OUTPUT_DIR, exist_ok=True)
    except OSError as e:
        logger.error(f'Error creating directory: {e}')

    logger.info(f'Created output directory: {HMMER_OUTPUT_DIR}')
    
    # setup the database and get some pairs to run
    tmpdir_database, db_path = PairPro.hmmer.create_accession_table()

    jaccard_threshold = 0.5
    chunk_size = 3
    
    PairPro.hmmer.process_pairs_table(db_path, chunk_size, HMMER_OUTPUT_DIR, jaccard_threshold)

    # creating model
    trained_model = model_dev(dbpath)[2]