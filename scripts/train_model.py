"""
Main script wrapper.
This is pseudocode for our parent script.

Directory: scipts

Input:

Output:

Note:
- We have two wrapping functions, but
these can be developed as individual scripts.
- Talk to Ryan about db_name and con
-

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
else:
    LOGLEVEL = 'DEBUG' # change to INFO for production
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'

##################
# Aux. functions #
##################
def combine_balanced_dfs(balanced_dfs, strategy='intersection'):
    """
    Combines multiple balanced dataframes based on a chosen strategy.

    Args:
        balanced_dfs (list): List of balanced dataframes.
        strategy (str): The combination strategy to use. Options are 'intersection' and 'union'.

    Returns:
        DataFrame: A combined dataframe.
    """

    if strategy == 'intersection':
        # Use merge to get the intersection of all dataframes
        combined_df = balanced_dfs[0]
        for df in balanced_dfs[1:]:
            combined_df = combined_df.merge(df, how='inner')
    
    elif strategy == 'union':
        # Use concat to get the union of all dataframes
        combined_df = pd.concat(balanced_dfs, axis=0).drop_duplicates().reset_index(drop=True)
    
    else:
        raise ValueError(f"Unknown combination strategy: {strategy}")

    return combined_df


def balance_dataframe(df, target_columns, strategy='undersample'):
    """
    Balances a dataframe based on the target column(s) and the chosen strategy.

    Args:
        df (DataFrame): The dataframe to balance.
        target_columns (list): List of target column names.
        strategy (str): The sampling strategy to use. Options are 'undersample', 'oversample', 'smote', and 'none'.

    Returns:
        list: A list of balanced dataframes for each target column.
    """
    # Ensure target_columns is a list, even if it's a single column.
    if not isinstance(target_columns, list):
        target_columns = [target_columns]

    balanced_dfs = []

    for target in target_columns:
        # separate the majority and minority classes
        majority_class = df[df[target] == True]
        minority_class = df[df[target] == False]

        if strategy == 'undersample':
            n_samples = len(minority_class)
            sampled_majority = resample(majority_class, n_samples=n_samples, replace=False)
            balanced_df = pd.concat([sampled_majority, minority_class])

        elif strategy == 'oversample':
            n_samples = len(majority_class)
            sampled_minority = resample(minority_class, n_samples=n_samples, replace=True)
            balanced_df = pd.concat([sampled_minority, majority_class])
        
        elif strategy == 'smote':
            raise NotImplementedError('SMOTE is not yet implemented')

        elif strategy == 'none':
            balanced_df = df

        else:
            raise ValueError(f"Unknown sampling strategy: {strategy}")

        balanced_dfs.append(balanced_df)
        # add logging statement to see the shape of balanced dataframe
        logger.debug(f'Dataframe shape after balancing for {target}: {balanced_df.shape}')

    return balanced_dfs # if combining mutliple target columns

################
# Main script #
################


@click.command()
@click.option('--hmmer', default=False, help='Whether to run HMMER or not')
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
def model_construction(hmmer, chunk_size, njobs, jaccard_threshold,
                       vector_size, structure, features):
    """
    Function to train a ML model to classify protein pairs
    """
    # build DB

    db_path = './tmp/pairpro.db'

    # connect to database
    con, _ = connect_db(TEST_DB_PATH)
    con, db_name = build_pairpro(con, db_path)
    db_path = f'./tmp/{db_name}.db'
    logger.info(f'Connected to database. Built pairpro table in {db_path}')

    ml_feature_list = []
    
    if hmmer:

        # press the HMM db
        pairpro.hmmer.hmmpress_hmms(HMM_PATH, PRESS_PATH)

        logger.info(f'Pressed HMM DB: {PRESS_PATH}')

        logger.info('Starting to run HMMER')

        # get all the proteins in pairs

        proteins_in_pair_count = con.execute(f"SELECT COUNT(*) FROM {db_name}.pairpro.proteins").fetchone()[0]
        # proteins_in_pair_count = con.execute(f"""SELECT COUNT(*) FROM (SELECT * FROM {db_name}.pairpro.proteins LIMIT 100) sub""").fetchone()[0]
        logger.debug(
            f"Total number of protein in pairs: {proteins_in_pair_count} in pipeline")

        proteins_in_pair = con.execute(
            f"SELECT pid, protein_seq FROM {db_name}.pairpro.proteins")
    
    
        # get number of hmms for evalue calc
        profiles = list(pyhmmer.plan7.HMMFile(HMM_PATH))
        n_hmms = len(profiles)
        del profiles
        logger.info(f"Number of HMMs: {n_hmms}")

        # run hmmsearch
        targets = pairpro.hmmer.prefetch_targets(PRESS_PATH)
        logger.debug(f"number of targets: {len(targets)}")
        wrapper = lambda chunk_index, pid_chunk: pairpro.hmmer.local_hmmer_wrapper(
            chunk_index, pid_chunk, press_path=PRESS_PATH, hmm_path=HMM_PATH, out_dir=HMMER_OUTPUT_DIR, cpu=njobs, prefetch=targets, e_value=1.e-5, scan=False, Z=n_hmms)

        complete = False
        chunk_index = 0
        total_processed = 0
        # use tqdm to track progress
        pbar = tqdm(total=proteins_in_pair_count)
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

            # update progress bar
            pbar.update(len(pid_chunk))
        pbar.close()
        
        logger.info('Starting to parse HMMER output')

        # setup the database and get some pairs to run
        con.execute("""
            CREATE OR REPLACE TABLE proteins_from_pairs AS
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
        con.execute("""CREATE OR REPLACE TEMP TABLE hmmer_results AS 
                    SELECT * FROM read_csv_auto('./data/protein_pairs/parsed_hmmer_output/*.csv', HEADER=TRUE)
                    WHERE functional IS NOT NULL AND score IS NOT NULL
                    """)
        con.execute(
            f"""ALTER TABLE {db_name}.pairpro.final ADD COLUMN hmmer_match BOOLEAN""")
        con.execute(f"""UPDATE {db_name}.pairpro.final AS f
        SET hmmer_match = hmmer.functional::BOOLEAN
        FROM hmmer_results AS hmmer
        WHERE 
            hmmer.meso_pid = f.meso_pid
            AND hmmer.thermo_pid = f.thermo_pid
            AND hmmer.functional IS NOT NULL
            AND hmmer.score IS NOT NULL;
        """)
        # delete rows from pairpro.final where corresponding hmmer_results have NaN functional
        # Delete rows from pairpro.final where hmmer_match is NULL
        con.execute(f"""DELETE FROM {db_name}.pairpro.final
        WHERE hmmer_match IS NULL;
        """)
        logger.info('Finished appending parsed HMMER output to table.')
        ml_feature_list.append('hmmer_match')

        df = con.execute(f"""SELECT pair_id, m_protein_seq, t_protein_seq, bit_score, local_gap_compressed_percent_id,
        scaled_local_query_percent_id, scaled_local_symmetric_percent_id,
        query_align_len, query_align_cov, subject_align_len, subject_align_cov,
        LENGTH(m_protein_seq) AS m_protein_len, LENGTH(t_protein_seq) AS t_protein_len, {', '.join(ml_feature_list)} FROM {db_name}.pairpro.final""").df()

        logger.debug(f"DataFrame shape after HMMER processing: {df.shape}")


    # structure component
    if structure:
        structure_df = con.execute(
            f"""SELECT pair_id, thermo_pid, thermo_pdb, meso_pid, meso_pdb FROM {db_name}.pairpro.final""").df()
        downloader = pairpro.structures.ProteinDownloader(pdb_dir=STRUCTURE_DIR)
        logger.info(
            f'Downloading structures. Output directory: {STRUCTURE_DIR}')
        downloader.download_structure(
            structure_df, 'meso_pdb', 'meso_pid')
        downloader.download_structure(
            structure_df, 'thermo_pdb', 'thermo_pid')
        logger.info('Finished downloading structures. Running FATCAT.')
        processor = pairpro.structures.FatcatProcessor(pdb_dir=STRUCTURE_DIR)
        processor.run_fatcat_dict_job(
            structure_df, output_file=f'{STRUCTURE_OUTPUT_DIR}output.csv')
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

        ml_feature_list.append('structure_match')

        df = con.execute(f"""SELECT pair_id, m_protein_seq, t_protein_seq, bit_score, local_gap_compressed_percent_id,
        scaled_local_query_percent_id, scaled_local_symmetric_percent_id,
        query_align_len, query_align_cov, subject_align_len, subject_align_cov,
        LENGTH(m_protein_seq) AS m_protein_len, LENGTH(t_protein_seq) AS t_protein_len, {', '.join(ml_feature_list)} FROM {db_name}.pairpro.final WHERE structure_match IS NOT NULL""").df()


        logger.debug(f"DataFrame shape after structure processing: {df.shape}")

    else:
        raise NotImplementedError('Currently, you cannot train a model without hmmer or structure')
        
    
    logger.debug(df.info(verbose=True))
    logger.debug(df.shape)
    logger.debug(df.head())
    logger.debug(df.keys())

    logger.info('Beginning to preprocess data for model training')

    # specify target (this is tenative for now. Will be updated later)
    if hmmer and structure:
        target = ['hmmer_match', 'structure_match']
    elif hmmer:
        target = 'hmmer_match'
    elif structure:
        target = 'structure_match'
    else:
        raise NotImplementedError('Currently, you cannot train a model without hmmer or structure')
    
    # balance the dataframe
    balanced_dataframes = balance_dataframe(df, target_columns=target, strategy='undersample')

    logger.debug(f"Number of balanced dataframes: {len(balanced_dataframes)}")
    logger.debug(f"DataFrame shape before balancing: {df.shape}")

    # combine the balanced dataframes
    df = combine_balanced_dfs(balanced_dataframes, strategy='union')

    logger.debug(f"DataFrame shape after balancing: {df.shape}")

    # you can use ifeature omega by enternig feature_list as feature
    if 'structure_match' in target:
        accuracy_score, model = train_val_wrapper(df, target, True, features)
        logger.info(f'Accuracy score: {accuracy_score}')

        joblib.dump(model, f'{MODEL_PATH}trained_model.pkl')
        logger.debug(f'model training data is {df.head()}')
        logger.info(f'Model saved to {MODEL_PATH}')
        con.close()

    else:
        accuracy_score, model = train_val_wrapper(df, target, False, features)
        logger.info(f'Accuracy score: {accuracy_score}')

        joblib.dump(model, f'{MODEL_PATH}trained_model.pkl')
        logger.debug(f'model training data is {df.head()}')
        logger.info(f'Model saved to {MODEL_PATH}')
        con.close()



        
if __name__ == "__main__":
    # Initialize logger
    logger = pairpro.utils.start_logger_if_necessary(
        LOGNAME, LOGFILE, LOGLEVEL, filemode='w')
    
    # Nested loggers
    hmmer_logger = logging.getLogger('pairpro.hmmer')
    hmmer_logger.setLevel(getattr(logging, LOGLEVEL))
    hmmer_logger.addHandler(logging.FileHandler(LOGFILE))

    structure_logger = logging.getLogger('pairpro.structures')
    structure_logger.setLevel(getattr(logging, LOGLEVEL))
    structure_logger.addHandler(logging.FileHandler(LOGFILE))


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
