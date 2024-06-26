"""
Python script to train the model on OMA data (migration of data)

Target: OMA label
Features: BLAST, HMMER, Structure, iFeatureOmega

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
import pairpro.utils as pp_utils
# hmmer
import pairpro.hmmer as pp_hmmer
# structure
import pairpro.structures as pp_structures
# ML
from pairpro.train_val_wrapper import train_val_wrapper


####################
### PATHS & VARS ###
####################
# db Paths
TEST_DB_PATH = '230830-truefalse_1k.csv' 

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

# Logan edit of combined dataframe function (need to change function call name in main script):

def balance_data(dataframe, target_columns):
    """
    Resamples the dataframe to evenly distribute labels

    Args:
        dataframe (pandas dataframe): training dataframe
        label_columns (list): list of columns to sample from

    Returns:
        pandas dataframe: New DF with evenly sampled labels
    """
    # Ensure target_columns is a list, even if it's a single column.
    if not isinstance(target_columns, list):
        target_columns = list(target_columns)

    for target in target_columns:
        # separate the majority and minority classes
        majority_class = dataframe[dataframe[target] == dataframe[target].value_counts().idxmax()]
        minority_class = dataframe[dataframe[target] == dataframe[target].value_counts().idxmin()]

        #create new dataframe with len(minority_class)
        n_samples = len(minority_class)
        undersampled_majority = resample(majority_class, n_samples=n_samples, replace=False)

        # Combine the undersampled majority class with the minority class
        dataframe = pd.concat([undersampled_majority, minority_class])
        logger.debug(f'DF length reduced to {dataframe.shape}')
        logger.debug(f'{target} value counts: {dataframe[target].value_counts()}')
        
    return dataframe

################
# Main script #
################


@click.command()
@click.option('--blast', default=False, help='Whether to run BLAST or not')
@click.option('--hmmer', default=False, help='Whether to run HMMER or not')
@click.option('--chunk_size', default=1,
              help='Number of sequences to process in each chunk. Size of 1 means 2048 sequences per chunk, 2 means 4096, etc. as it is vectorized.')
@click.option('--njobs', default=4,
              help='Number of parallel processes to use for HMMER')
@click.option('--jaccard_threshold', default=0.7,
              help='Jaccard threshold for filtering protein pairs')
@click.option('--features', default=False,
              help='List of features to use for the model')
@click.option('--structure', default=False,
              help='Whether to use structure or not')
@click.option('--model_name', default='base',
              help='Name of the model')
def model_construction(blast, hmmer, chunk_size, njobs, jaccard_threshold,
                       structure, features, model_name):
    """
    Function to train a ML model to classify protein pairs
    """

    ##### database construction #####

    con = ddb.connect('./tmp/OMA_sample.db', read_only=False) # create a database. Has to be read_only=False

    # create a main table
    con.execute("""CREATE OR REPLACE TABLE OMA_1k AS SELECT query, subject, query_id, subject_id, pair_id, 
                local_gap_compressed_percent_id, scaled_local_query_percent_id, scaled_local_symmetric_percent_id, query_align_len, 
                query_align_cov, subject_align_len, subject_align_cov, bit_score, Pair, jaccard FROM read_csv_auto('230830-truefalse_1k.csv')""")
    con.commit() # commit the changes. Otherwise, the table will not be created.
    # create a table for proteins in pairs
    con.execute("""CREATE OR REPLACE TABLE proteins AS 
        (
            SELECT DISTINCT pid, protein_seq
            FROM 
            (
                SELECT query_id AS pid, query as protein_seq
                FROM OMA_1k
                UNION ALL
                SELECT subject_id AS pid, subject as protein_seq
                FROM OMA_1k
            )   
        );""")


    # # build DB

    # db_path = './tmp/pairpro.db'

    # # connect to database
    # con, _ = connect_db(TEST_DB_PATH)
    # con, db_name = build_pairpro(con, db_path)
    # db_path = f'./tmp/{db_name}.db'
    # logger.info(f'Connected to database. Built pairpro table in {db_path}')

    ml_feature_list = [] # list of features to use for ML

    ##### feature construction #####

    if blast:
        # blast component
        logger.info('Starting to run BLAST')
        pass
        # raise NotImplementedError("BLAST is not implemented yet") ## uncomment this when switching away from OMA_1k
    else:
        pass
    
    if hmmer: #pids were strings in the original code

        # press the HMM db
        pp_hmmer.hmmpress_hmms(HMM_PATH, PRESS_PATH)

        logger.info(f'Pressed HMM DB: {PRESS_PATH}')

        logger.info('Starting to run HMMER')

        # get all the proteins in pairs

        proteins_in_pair_count = con.execute(f"SELECT COUNT(*) FROM proteins").fetchone()[0]
        # proteins_in_pair_count = con.execute(f"""SELECT COUNT(*) FROM (SELECT * FROM {db_name}.pairpro.proteins LIMIT 100) sub""").fetchone()[0]
        logger.debug(
            f"Total number of protein in pairs: {proteins_in_pair_count} in pipeline")

        proteins_in_pair = con.execute(
            f"SELECT pid, protein_seq FROM proteins")
    
    
        # get number of hmms for evalue calc
        profiles = list(pyhmmer.plan7.HMMFile(HMM_PATH))
        n_hmms = len(profiles)
        del profiles
        logger.info(f"Number of HMMs: {n_hmms}")

        # run hmmsearch
        targets = pp_hmmer.prefetch_targets(PRESS_PATH)
        logger.debug(f"number of targets: {len(targets)}")
        wrapper = lambda chunk_index, pid_chunk: pp_hmmer.local_hmmer_wrapper(
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
        pp_hmmer.process_pairs_table(
            con,
            chunk_size,
            PARSE_HMMER_OUTPUT_DIR,
            jaccard_threshold)
        

        logger.info('Finished parsing HMMER output.')

        # checking if the parsed output is appended to table
        con.execute("""CREATE OR REPLACE TEMP TABLE hmmer_results AS 
                    SELECT * FROM read_csv_auto('./data/protein_pairs/parsed_hmmer_output/*.csv', HEADER=TRUE)
                    WHERE functional IS NOT NULL AND score IS NOT NULL
                    """)
        con.execute(
            f"""ALTER TABLE OMA_1k ADD COLUMN hmmer_match DOUBLE""")
        con.execute(f"""UPDATE OMA_1k AS f
        SET hmmer_match = hmmer.score::DOUBLE
        FROM hmmer_results AS hmmer
        WHERE 
            hmmer.query_id = f.query_id
            AND hmmer.subject_id = f.subject_id
            AND hmmer.functional IS NOT NULL
            AND hmmer.score IS NOT NULL;
        """)
        # delete rows from pairpro.final where corresponding hmmer_results have NaN functional
        # Delete rows from pairpro.final where hmmer_match is NULL
        con.execute(f"""DELETE FROM OMA_1k
        WHERE hmmer_match IS NULL;
        """)
        logger.info('Finished appending parsed HMMER output to table.')
        ml_feature_list.append('hmmer_match')

        df = con.execute(f"""SELECT pair_id, query, subject, bit_score, local_gap_compressed_percent_id,
        scaled_local_query_percent_id, scaled_local_symmetric_percent_id,
        query_align_len, query_align_cov, subject_align_len, subject_align_cov,
        LENGTH(query) AS query_len, LENGTH(subject) AS subject_len, Pair, {', '.join(ml_feature_list)} FROM OMA_1k""").df()

        logger.debug(f"DataFrame shape after HMMER processing: {df.shape}")
    else:
        pass


    # structure component ###leave this for now. Check with Ryan on OMA database
    if structure:
        structure_df = con.execute(
            f"""SELECT pair_id, thermo_pid, thermo_pdb, meso_pid, meso_pdb FROM OMA_1k""").df()
        downloader = pp_structures.ProteinDownloader(pdb_dir=STRUCTURE_DIR)
        logger.info(
            f'Downloading structures. Output directory: {STRUCTURE_DIR}')
        downloader.download_structure(
            structure_df, 'meso_pdb', 'meso_pid')
        downloader.download_structure(
            structure_df, 'thermo_pdb', 'thermo_pid')
        logger.info('Finished downloading structures. Running FATCAT.')
        processor = pp_structures.FatcatProcessor(pdb_dir=STRUCTURE_DIR)
        processor.run_fatcat_dict_job(
            structure_df, output_file=f'{STRUCTURE_OUTPUT_DIR}output.csv')
        logger.info('Finished running FATCAT.')

        con.execute("""CREATE OR REPLACE TEMP TABLE structure_results AS SELECT * FROM read_csv_auto('./data/protein_pairs/structures/*.csv', HEADER=TRUE)""")
        con.execute(
            f"""ALTER TABLE OMA_1k ADD COLUMN structure_match BOOLEAN""")
        con.execute(f"""UPDATE OMA_1k AS f
        SET structure_match = structure.p_value::BOOLEAN
        FROM structure_results AS structure
        WHERE structure.pair_id = f.pair_id
        """)
        logger.info('Finished appending structure output to table.')

        ml_feature_list.append('structure_match')

        df = con.execute(f"""SELECT pair_id, query, subjecy, bit_score, local_gap_compressed_percent_id,
        scaled_local_query_percent_id, scaled_local_symmetric_percent_id,
        query_align_len, query_align_cov, subject_align_len, subject_align_cov,
        LENGTH(query) AS query_len, LENGTH(subject) AS subject_len, Pair, {', '.join(ml_feature_list)} FROM OMA_1k WHERE structure_match IS NOT NULL""").df()


        logger.debug(f"DataFrame shape after structure processing: {df.shape}")

    else:
        pass

    if not hmmer and not structure:
        logger.info('No features selected. Exiting.')
        raise NotImplementedError('Currently, you cannot train a model without hmmer or structure')
    else:
        pass
        
    
    logger.debug(df.info(verbose=True))
    logger.debug(df.shape)
    logger.debug(df.head())
    logger.debug(df.keys())

    logger.info('Beginning to preprocess data for model training')

    # # specify target (this is tenative for now. Will be updated later)
    # if hmmer and structure:
    #     target = ['hmmer_match', 'structure_match']
    # elif hmmer:
    #     target = 'hmmer_match'
    # elif structure:
    #     target = 'structure_match'
    # else:
    #     raise NotImplementedError('Currently, you cannot train a model without hmmer or structure')

    target = ['Pair']
    
    # balance the dataframe (Logan version)
    df = balance_data(df, target_columns=target)
    logger.debug(f"DataFrame shape after balancing: {df.shape}")


    # # you can use ifeature omega by enternig feature_list as feature
    # if 'structure_match' in target:
    #     accuracy_score, model = train_val_wrapper(df, target, blast, hmmer, structure, features)
    #     logger.info(f'Accuracy score: {accuracy_score}')

    #     joblib.dump(model, f'{MODEL_PATH}{model_name}.pkl')
    #     logger.debug(f'model training data is {df.head()}')
    #     logger.info(f'Model saved to {MODEL_PATH}')
    #     con.close()
    # else:
    accuracy_score, model = train_val_wrapper(df, target, blast, hmmer, structure, features)
    logger.info(f'Accuracy score: {accuracy_score}')

    joblib.dump(model, f'{MODEL_PATH}{model_name}.pkl')
    logger.debug(f'model training data is {df.head()}')
    logger.info(f'Model saved to {MODEL_PATH}')
    con.close()



        
if __name__ == "__main__":
    # Initialize logger
    logger = pp_utils.start_logger_if_necessary(
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