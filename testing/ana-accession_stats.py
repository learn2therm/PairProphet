"""
This scripts aims to generate statistics and plots for the accession data generated by HMMER.
In particular, it will do this by varing two main parameters: 1) e-value and 2) Jaccard similarity threshold.
These two parameters will be varied in order to determine the optimal values for each (if such a value exists).
And morevoer, have a stronger understanding of the distribution of the data as it goes through the pipeline 
especially machine learning.

Environment Variables to set:
    LOGLEVEL: logging level for the script (default: INFO, options: DEBUG, INFO, WARNING, ERROR, CRITICAL)

Note:
    This is a script meant to run at base level of the repository. e.g. 'python testing/ana-accession_stats.py'
    It also assumes that 'train_model.py' (and associated scripts) have been run and the data is in the correct location.
    In general, I assume the 50k dataset is used for this analysis. So if you want to use a different dataset, you will need to change the paths
    and db names accordingly.

Future Work:
    - [ ] make paths options with default values
    - [ ] LOGLEVEL should be in the CLI (and not env var)
    - [ ] think about a better data structure for evalue and jaccard threshold values  
    - [ ] use os more for paths stuff specifically for assembling paths
        os.path.join(path1, path2, path3, ...)
"""
# system dependencies
import os
import logging

# library dependencies
import click
import duckdb as ddb

import pandas as pd
import pyhmmer
from tqdm import tqdm


# local dependencies
from pairpro.preprocessing import connect_db, build_pairpro
from pairpro.user_blast import make_blast_df

import pairpro.hmmer as pp_hmmer
import pairpro.utils as pp_utils

####################
### PATHS & VARS ###
####################

# Paths
HMM_PATH = './data/pfam/Pfam-A.hmm'
DB_PATH = './tmp/pairpro.db'
PRESS_PATH = './data/pfam/pfam'

ANALYSIS_OUTPUT_PATH = './data/analysis/'
HMMER_OUTPUT_DIR = './data/analysis/hmmer/'
PARSE_HMMER_OUTPUT_DIR = './data/analysis/hmmer_parsed/'


# venv variables
if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
    LOGLEVEL = getattr(logging, LOGLEVEL)
else:
    LOGLEVEL = logging.INFO
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'

####################



####################
###    SCRIPT    ###
####################


@click.command()
@click.option('--chunk_size', default=1,
              help='Number of sequences to process in each chunk')
@click.option('--njobs', default=4,
              help='Number of parallel processes to use for HMMER')
@click.option('--evalue', default=1.e-10,
                help='E-value for HMMER')
@click.option('--jaccard_threshold', default=0.7,
              help='Jaccard threshold for filtering protein pairs')
@click.option('--vector_size', default=1,
              help='Size of the vector for the dataframe chunking') 
def analysis_script(chunk_size, njobs, evalue, jaccard_threshold, vector_size, **kwargs):
    """
    This a wrapper function that will run the analysis for the accession data.
    """
    # connect to database
    con = ddb.connect(DB_PATH, read_only=False)

    # get number of proteins in pairs
    proteins_in_pair_count = con.execute(f"SELECT COUNT(*) FROM pairpro.pairpro.proteins").fetchone()[0]
    logger.debug(
        f"Total number of protein in pairs: {proteins_in_pair_count} in pipeline")

    # get proteins in pairs
    proteins_in_pair = con.execute(
        f"SELECT pid, protein_seq FROM pairpro.pairpro.proteins")
    
    # create empty lists to store statistics
    evalue_values = []
    jaccard_threshold_values = []
    mean_acc_length_values = []

    for evalue_value in evalue_values_to_test:
        for jaccard_threshold_value in jaccard_threshold_values_to_test:

            ### Run HMMER ###
            # get number of hmms for evalue calc
            profiles = list(pyhmmer.plan7.HMMFile(HMM_PATH))
            n_hmms = len(profiles)
            del profiles
            logger.info(f"Number of HMMs: {n_hmms}")

            # run hmmsearch
            targets = pp_hmmer.prefetch_targets(PRESS_PATH)
            logger.info(f"Number of targets: {len(targets)}")
            wrapper = lambda chunk_index, pid_chunk: pp_hmmer.local_hmmer_wrapper(
                chunk_index, pid_chunk, press_path=PRESS_PATH, hmm_path=HMM_PATH, out_dir=HMMER_OUTPUT_DIR, cpu=njobs, prefetch=targets, e_value=evalue, scan=False, Z=n_hmms)
            
            # the hmmsearch loop
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

            ### Read HMMER output ###
            logger.info('Reading HMMER output...')
            df_lists = [] # initialize list to store dataframes
            for chunk_index in range(chunk_index):
                output_file_path = f'{HMMER_OUTPUT_DIR}{chunk_index}_output.csv'
                df = pd.read_csv(output_file_path)
                logger.debug(f"Loaded chunk {chunk_index} with size {len(df)}")
                df_lists.append(df)
                

            logger.debug(f"df_list: {df_lists}")
            bfd = pd.concat(df_lists, ignore_index=True)
        
            # calculate stats on the data (e.g. mean length of accessions)
            mean_accession_length = bfd['accession_id'].str.split(';').str.len().mean()
            
            # .apply(lambda x: len(x.split(';'))).mean()
            logger.debug(f"Mean length of accessions: {mean_accession_length}")            

            # store the stats for each combination
            evalue_values.append(evalue_value)
            jaccard_threshold_values.append(jaccard_threshold_value)
            mean_acc_length_values.append(mean_accession_length)
    
    logger.info('Finished running HMMER. Added stats')

    # store the stats in a dataframe
    results_df = pd.DataFrame({
        'e-value': evalue_values,
        'jaccard_threshold': jaccard_threshold_values,
        'mean_acc_length': mean_acc_length_values
    })

    # save the result in a dataframe for later use
    results_df.to_csv(f'{ANALYSIS_OUTPUT_PATH}statistics_results.csv', index=False)

    # parse the output
    logger.info('Parsing HMMER output...')

    # # setup the database and get some pairs to run
    # con.execute("""
    #     CREATE TABLE proteins_from_pairs4 AS
    #     SELECT query_id AS pid, accession_id AS accession
    #     FROM read_csv_auto('./data/analysis/hmmer/*.csv', HEADER=TRUE)
    # """)
    # con.commit()

    logger.info('creating pair table')
    pp_hmmer.process_pairs_table_ana(
        con,
        'pairpro',
        vector_size,
        PARSE_HMMER_OUTPUT_DIR,
        jaccard_threshold)

    logger.info('Finished parsing HMMER output.')

    


if __name__ == "__main__":
    # Initialize logger
    logger = pp_utils.start_logger_if_necessary(
        LOGNAME, LOGFILE, LOGLEVEL, filemode='w')
    logger.info(f"Running {__file__}")

    # create analysis directory
    try:
        os.makedirs(ANALYSIS_OUTPUT_PATH, exist_ok=True)
    except OSError as e:
        logger.error(f'Error creating directory: {e}')
    
    # create HMMER output directory
    try:
        os.makedirs(HMMER_OUTPUT_DIR, exist_ok=True)
    except OSError as e:
        logger.error(f'Error creating directory: {e}')

    # create parsed HMMER output directory
    try:
        os.makedirs(PARSE_HMMER_OUTPUT_DIR, exist_ok=True)
    except OSError as e:
        logger.error(f'Error creating directory: {e}')


    ## assume that the data is already in the correct location
    ## the hmm's have been also pressed

    # define the range of e-value and Jaccard threshold values to test
    evalue_values_to_test = [1e-10, 1e-5, 1e-3, 1e-1]
    jaccard_threshold_values_to_test = [0.4, 0.5, 0.7, 0.9]


    # run the analysis
    analysis_script()

    
