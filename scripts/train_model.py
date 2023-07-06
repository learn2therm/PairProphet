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
This script will take a long time to run.
"""
# system dependencies
import sys
import logging
import os

# library dependencies
import click
import duckdb as ddb
import pandas as pd
import pyhmmer
import joblib
from joblib import Parallel, delayed
from sklearn.utils import resample
from tqdm import tqdm

# local dependencies

# machine learning
# from pairpro.evaluate_model import evaluate_model
# from pairpro.train_val_featuregen import create_new_dataframe
from pairpro.train_val_wrapper import train_val_wrapper
# from pairpro.train_val_input_cleaning import columns_to_keep

# build DB
from pairpro.preprocessing import connect_db, build_pairpro
from pairpro.user_blast import make_blast_df

# hmmer
import pairpro.hmmer
import pairpro.utils


# structure
import pairpro.structures


# db Paths
TEST_DB_PATH = './data/l2t_50k.db'  # l2t_50k.db

# HMMER Paths
HMM_PATH = './data/pfam/Pfam-A.hmm'  # ./Pfam-A.hmm
PRESS_PATH = './data/pfam/pfam'
HMMER_OUTPUT_DIR = './data/protein_pairs/'
PARSE_HMMER_OUTPUT_DIR = './data/protein_pairs/parsed_hmmer_output/'
WORKER_WAKE_UP_TIME = 25  # this is to ensure that if a worker that is about to be shut down due to previous task completetion doesn't actually start running

# Structure Paths
STRUCTURE_DIR = './data/structures/'
STRUCTURE_OUTPUT_DIR = './data/protein_pairs/structures/'

# ML Paths
MODEL_PATH = './data/models/'

# get environmental variables
if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
    LOGLEVEL = getattr(logging, LOGLEVEL)
else:
    LOGLEVEL = logging.INFO
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'


@click.command()
@click.option('--chunk_size', default=1,
              help='Number of sequences to process in each chunk')
@click.option('--njobs', default=4,
              help='Number of parallel processes to use for HMMER')
@click.option('--jaccard_threshold', default=0.7,
              help='Jaccard threshold for filtering protein pairs')
@click.option('--vector_size', default=1,
              help='Size of the vector for the dataframe chunking')
@click.option('--features', default=False,
              help='List of features to use for the model')
@click.option('--structure', default=False,
              help='Whether to use structure or not')
def model_construction(chunk_size, njobs, jaccard_threshold,
                       vector_size, structure, features):
    """_summary_
    """
    # press the HMM db
    pairpro.hmmer.hmmpress_hmms(HMM_PATH, PRESS_PATH)

    logger.info(f'Pressed HMM DB: {PRESS_PATH}')

    db_path = './tmp/pairpro.db'

    # connect to database
    con, _ = connect_db(TEST_DB_PATH)
    con, db_name = build_pairpro(con, db_path)
    db_path = f'./tmp/{db_name}.db'
    logger.info(f'Connected to database. Built pairpro table in {db_path}')

    logger.info('Starting to run HMMER')

    # get all the proteins in pairs

    proteins_in_pair = con.execute(
        f"SELECT pid, protein_seq FROM {db_name}.pairpro.proteins") #take out df() later
    # logger.debug(
    #     f"Total number of protein in pairs: {len(proteins_in_pair)} in pipeline")
    
    # get number of hmms for evalue calc
    profiles = list(pyhmmer.plan7.HMMFile(HMM_PATH))
    n_hmms = len(profiles)
    del profiles
    logger.info(f"Number of HMMs: {n_hmms}")

    # run hmmsearch
    targets = pairpro.hmmer.prefetch_targets(PRESS_PATH)
    logger.debug(f"number of targets: {len(targets)}")
    wrapper = lambda chunk_index, pid_chunk: pairpro.hmmer.local_hmmer_wrapper(
        chunk_index, pid_chunk, press_path=PRESS_PATH, hmm_path=HMM_PATH, out_dir=HMMER_OUTPUT_DIR, cpu=njobs, prefetch=targets, e_value=1.e-10, scan=False, Z=n_hmms)

    complete = False
    chunk_index = 0
    total_processed = 0
    while not complete:
        pid_chunk = proteins_in_pair.fetch_df_chunk(vectors_per_chunk=chunk_size)
        logger.info(f"Loaded chunk of size {len(pid_chunk)}")
        if len(pid_chunk) == 0:
            complete = True
            break
        wrapper(chunk_index, pid_chunk)
        logger.info(f"Ran chunk, validating results")

        df = pd.read_csv(f'{HMMER_OUTPUT_DIR}/{chunk_index}_output.csv')
        assert set(list(pid_chunk['pid'].values)) == set(list(df['query_id'].values)), "Not all query ids are in the output file"

        logger.info(f"Completed chunk {chunk_index} with size {len(pid_chunk)}")
        total_processed += len(pid_chunk)
        chunk_index += 1

    
    logger.info('Starting to parse HMMER output')

    # setup the database and get some pairs to run
    con.execute("""
        CREATE TABLE proteins_from_pairs AS
        SELECT query_id AS pid, accession_id AS accession
        FROM read_csv_auto('./data/protein_pairs/*.csv', HEADER=TRUE)
    """)
    con.commit()

    logger.info(
        f'Created parse HMMER output directory: {PARSE_HMMER_OUTPUT_DIR}. Running parse HMMER algorithm.')
    pairpro.hmmer.process_pairs_table(
        con,
        db_name,
        vector_size,
        PARSE_HMMER_OUTPUT_DIR,
        jaccard_threshold)

    logger.info('Finished parsing HMMER output.')

    # checking if the parsed output is appended to table
    con.execute("""CREATE TABLE hmmer_results AS SELECT * FROM read_csv_auto('./data/protein_pairs/parsed_hmmer_output/*.csv', HEADER=TRUE)""")
    con.execute(
        f"""ALTER TABLE {db_name}.pairpro.final ADD COLUMN hmmer_match BOOLEAN""")
    con.execute(f"""UPDATE {db_name}.pairpro.final AS f
    SET hmmer_match = hmmer.functional::BOOLEAN
    FROM hmmer_results AS hmmer
    WHERE hmmer.meso_pid = f.meso_pid
    AND hmmer.thermo_pid = f.thermo_pid
    """)
    logger.info('Finished appending parsed HMMER output to table.')

    # structure component
    if structure:
        structure_df = con.execute(
            f"""SELECT pair_id, thermo_pid, thermo_pdb, meso_pid, meso_pdb FROM {db_name}.pairpro.final""").df()
        logger.info(
            f'Downloading structures. Output directory: {STRUCTURE_DIR}')
        pairpro.structures.download_structure(
            structure_df, 'meso_pdb', 'meso_pid', STRUCTURE_DIR)
        pairpro.structures.download_structure(
            structure_df, 'thermo_pdb', 'thermo_pid', STRUCTURE_DIR)
        logger.info('Finished downloading structures. Running FATCAT.')
        pairpro.structures.run_fatcat_dict_job(
            structure_df, STRUCTURE_DIR, f'{STRUCTURE_OUTPUT_DIR}/output.csv')
        logger.info('Finished running FATCAT.')

        con.execute("""CREATE OR REPLACE TEMP TABLE structure_results AS SELECT * FROM read_csv_auto('./data/protein_pairs/structures/*.csv', HEADER=TRUE)""")
        con.execute(
            f"""ALTER TABLE {db_name}.pairpro.final ADD COLUMN structure_match BOOLEAN""")
        con.execute(f"""UPDATE {db_name}.pairpro.final AS f
        SET structure_match = structure.p_value::BOOLEAN
        FROM structure_results AS structure
        WHERE structure.pair_id = f.pair_id
        """)
        logger.info('Finished appending structure output to table.')

        df = con.execute(f"""SELECT pair_id, m_protein_seq, t_protein_seq, bit_score, local_gap_compressed_percent_id,
        scaled_local_query_percent_id, scaled_local_symmetric_percent_id,
        query_align_len, query_align_cov, subject_align_len, subject_align_cov,
        LENGTH(m_protein_seq) AS m_protein_len, LENGTH(t_protein_seq) AS t_protein_len, hmmer_match, structure_match FROM {db_name}.pairpro.final WHERE structure_match IS NOT NULL""").df()

    else:
        logger.info('Skipping structure component.')
        df = con.execute(f"""SELECT pair_id, m_protein_seq, t_protein_seq, bit_score, local_gap_compressed_percent_id,
        scaled_local_query_percent_id, scaled_local_symmetric_percent_id,
        query_align_len, query_align_cov, subject_align_len, subject_align_cov,
        LENGTH(m_protein_seq) AS m_protein_len, LENGTH(t_protein_seq) AS t_protein_len, hmmer_match FROM {db_name}.pairpro.final""").df()

    print(df.head())
    print(df.shape)

    logger.info('Beginning to preprocess data for model training')

    # Separate the majority and minority classes
    majority_class = df[df['hmmer_match'] == True]
    minority_class = df[df['hmmer_match'] == False]

    # Undersample the majority class to match the number of minority class
    # samples
    n_samples = len(minority_class)
    undersampled_majority = resample(
        majority_class,
        n_samples=n_samples,
        replace=False)

    # Combine the undersampled majority class with the minority class
    df = pd.concat([undersampled_majority, minority_class])

    # specify hmmer target
    if structure:
        target = ['hmmer_match', 'structure_match']
    else:
        target = 'hmmer_match'

    # you can use ifeature omega by enternig feature_list as feature
    accuracy_score, model = train_val_wrapper(df, target, structure, features)
    logger.info(f'Accuracy score: {accuracy_score}')

    joblib.dump(model, f'{MODEL_PATH}trained_model.pkl')
    logger.debug(f'model training data is {df.head()}')
    logger.info(f'Model saved to {MODEL_PATH}')


if __name__ == "__main__":
    # Initialize logger
    logger = pairpro.utils.start_logger_if_necessary(
        LOGNAME, LOGFILE, LOGLEVEL, filemode='w')
    logger.info(f"Running {__file__}")

    # create pfam HMM directory (this was before HMM download script)
    try:
        os.makedirs('./data/pfam', exist_ok=True)
    except OSError as e:
        logger.error(f'Error creating directory: {e}')

    # prepare output file
    try:
        os.makedirs(HMMER_OUTPUT_DIR, exist_ok=True)
    except OSError as e:
        logger.error(f'Error creating directory: {e}')

    logger.info(f'Created HMMER output directory: {HMMER_OUTPUT_DIR}')

    # prepare output file
    try:
        os.makedirs(PARSE_HMMER_OUTPUT_DIR, exist_ok=True)
    except OSError as e:
        logger.error(f'Error creating directory: {e}')

    # prepare output file
    try:
        os.makedirs(STRUCTURE_DIR, exist_ok=True)
    except OSError as e:
        logger.error(f'Error creating directory: {e}')

    # prepare output file
    try:
        os.makedirs(STRUCTURE_OUTPUT_DIR, exist_ok=True)
    except OSError as e:
        logger.error(f'Error creating directory: {e}')

    # prepare output file
    try:
        os.makedirs(MODEL_PATH, exist_ok=True)
    except OSError as e:
        logger.error(f'Error creating directory: {e}')
    logger.info(f'Created model directory: {MODEL_PATH}')

    # creating model
    model_construction()
