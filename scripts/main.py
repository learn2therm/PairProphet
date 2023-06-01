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
TODO:
    - [x] Check HMMER parsing logic
"""
# system dependencies
import sys
import logging
import os

# library dependencies
import duckdb as ddb
import pandas as pd
import joblib
from joblib import Parallel, delayed
from sklearn.utils import resample
from tqdm import tqdm

# local dependencies

##machine learning
from PairPro.evaluate_model import evaluate_model
from PairPro.train_val_wrapper import train_val_wrapper
from PairPro.train_val_input_cleaning import columns_to_keep

##build DB
from PairPro.preprocessing import connect_db, build_pairpro
from PairPro.user_blast import make_blast_df

##hmmer
import PairPro.hmmer
import PairPro.utils


##structure
# from PairPro.structure import download_structures, run_fatcat


## db Paths
TEST_DB_PATH = '/Users/humoodalanzi/pfam/l2t_50k.db'

## HMMER Paths
HMM_PATH = './data/pfam/Pfam-A.hmm'  # ./Pfam-A.hmm
PRESS_PATH = './data/pfam/pfam'
HMMER_OUTPUT_DIR = './data/protein_pairs/'
PARSE_HMMER_OUTPUT_DIR = './data/protein_pairs/parsed_hmmer_output/'
WORKER_WAKE_UP_TIME = 25 # this is to ensure that if a worker that is about to be shut down due to previous task completetion doesn't actually start running

## ML Paths
MODEL_PATH = './data/models/'

## get environmental variables
if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
    LOGLEVEL = getattr(logging, LOGLEVEL)
else:
    LOGLEVEL = logging.INFO
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'


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

    db_path = './tmp/pairpro.db'

    #connect to database
    con, _ = connect_db(TEST_DB_PATH)
    con = build_pairpro(con, db_path)
    logger.info(f'Connected to database. Built pairpro table in {db_path}')

    logger.info('Starting to run HMMER')

    # Set up parallel processing and parsing
    chunk_size = 2500 # Number of sequences to process in each chunk 
    njobs = 4  # Number of parallel processes to use

    logger.info('Parallel processing parameters obtained')

    # prepare output file
    try:
        os.makedirs(HMMER_OUTPUT_DIR, exist_ok=True)
    except OSError as e:
        logger.error(f'Error creating directory: {e}')

    logger.info(f'Created HMMER output directory: {HMMER_OUTPUT_DIR}')

    proteins_in_pair_pids = con.execute("SELECT pid FROM pairpro.pairpro.proteins").df()
    logger.debug(f"Total number of protein in pairs: {len(proteins_in_pair_pids)} in pipeline")
    
    # chunking the PID so the worker function queries
    protein_pair_pid_chunks = [proteins_in_pair_pids[i:i + chunk_size]
                            for i in range(0, len(proteins_in_pair_pids), chunk_size)]
    
    con.close()
    
    # run hmmscan
    logger.info('Running pyhmmer in parallel on all chunks')

    with tqdm(total=len(protein_pair_pid_chunks)) as pbar:
        Parallel(
            n_jobs=njobs)(
            delayed(PairPro.hmmer.local_hmmer_wrapper)(
                chunk_index,
                db_path,
                protein_pair_pid_chunks,
                PRESS_PATH,
                HMMER_OUTPUT_DIR,
                None) for chunk_index,
            protein_pair_pid_chunks in enumerate(protein_pair_pid_chunks))
        pbar.update(1)
    logger.debug(f"number of protein chunks: {len(protein_pair_pid_chunks)}")

    logger.info('Finished running pyhmmer in parallel on all chunks')
    logger.info('Starting to parse HMMER output')
    
    # re-open connection (ask Ryan about this)
    con = ddb.connect(db_path)

    # setup the database and get some pairs to run
    con.execute("""
        CREATE TABLE proteins_from_pairs AS
        SELECT query_id AS pid, accession_id AS accession
        FROM read_csv_auto('./data/protein_pairs/*.csv', HEADER=TRUE)
    """)
    con.commit()

    # prepare output file
    try:
        os.makedirs(PARSE_HMMER_OUTPUT_DIR, exist_ok=True)
    except OSError as e:
        logger.error(f'Error creating directory: {e}')

    jaccard_threshold = 0.7
    chunk_size = 3
    
    logger.info(f'Created parse HMMER output directory: {PARSE_HMMER_OUTPUT_DIR}. Running parse HMMER algorithm.')
    PairPro.hmmer.process_pairs_table(con, chunk_size, PARSE_HMMER_OUTPUT_DIR, jaccard_threshold)

    logger.info('Finished parsing HMMER output.')

    # checking if the parsed output is appended to table
    con.execute("""CREATE TABLE hmmer_results AS SELECT * FROM read_csv_auto('./data/protein_pairs/parsed_hmmer_output/*.csv', HEADER=TRUE)""")
    con.execute("""ALTER TABLE pairpro.pairpro.final ADD COLUMN hmmer_match BOOLEAN""")
    con.execute("""UPDATE pairpro.pairpro.final AS f
    SET hmmer_match = hmmer.functional::BOOLEAN
    FROM hmmer_results AS hmmer
    WHERE hmmer.meso_pid = f.meso_pid
    AND hmmer.thermo_pid = f.thermo_pid
    """)
    logger.info('Finished appending parsed HMMER output to table.')

    df = con.execute("""SELECT bit_score, local_gap_compressed_percent_id, 
    scaled_local_query_percent_id, scaled_local_symmetric_percent_id, 
    query_align_len, query_align_cov, subject_align_len, subject_align_cov, 
    LENGTH(m_protein_seq) AS m_protein_len, LENGTH(t_protein_seq) AS t_protein_len, hmmer_match FROM pairpro.pairpro.final""").df()

    print(df.head())
    print(df.shape)

    logger.info('Beginning to preprocess data for model training')

    # Separate the majority and minority classes
    majority_class = df[df['hmmer_match'] == True]
    minority_class = df[df['hmmer_match'] == False]

    # Undersample the majority class to match the number of minority class samples
    n_samples = len(minority_class)
    undersampled_majority = resample(
        majority_class,
        n_samples=n_samples,
        replace=False)

    # Combine the undersampled majority class with the minority class
    df = pd.concat([undersampled_majority, minority_class])

    # creating model
    # prepare output file
    try:
        os.makedirs(MODEL_PATH, exist_ok=True)
    except OSError as e:
        logger.error(f'Error creating directory: {e}')
    logger.info(f'Created model directory: {MODEL_PATH}')

    accuracy_score = train_val_wrapper(df)[0]
    model = train_val_wrapper(df)[1]
    logger.info(f'Accuracy score: {accuracy_score}')

    joblib.dump(model, f'{MODEL_PATH}trained_model.pkl')
    logger.debug(f'model training data is {df.head()}')
    logger.info(f'Model saved to {MODEL_PATH}')
    