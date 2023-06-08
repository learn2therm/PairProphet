"""
This script is for user to interact with the pretrained model.

TODO:
- [X] Add click functionality
- [ ] Once structure component is ready, add structure component to the script
    - note: structure component has been given issues with Hyak, therefore, we did not include it in this iteration of the script.
            However, we run the structure component locally in a windows machine and it works, but it takes a long time to run.
- [ ] Add feature generation component to the script
"""

import sys
import logging
import os

# library dependencies
import click
import pandas as pd
import numpy as np
import sklearn
import joblib
from joblib import delayed, Parallel

# local dependencies
# machine learning
from pairpro.evaluate_model_wrapper import evaluate_model_wrapper

# need to understand how to import the trained model from main
# from pairpro.main import train_model

# build DB
from pairpro.user_blast import make_blast_df

# hmmer
import pairpro.hmmer
import pairpro.utils


# structure
# from pairpro.structures import download_structure, run_fatcat

# Paths
# ML Paths
MODEL_PATH = './data/models/trained_model.pkl'

# HMMER Paths
PRESS_PATH = './data/pfam/pfam'
HMMER_OUT_DIR = './data/user/hmmer_out/'
PARSED_HMMER_OUT_DIR = './data/user/parsed_hmmer_out/'

test_sequences = './data/50k_paired_seq.csv'


# get environmental variables
if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
    LOGLEVEL = getattr(logging, LOGLEVEL)
else:
    LOGLEVEL = logging.INFO
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'


@click.command()
@click.option('--test_sequences', default='./data/50k_paired_seq.csv',
              help='Path to csv with test sequences')
@click.option('--structure', default=False,
              help='Boolean; Run structure component')
@click.option('--features', default=False,
              help='Boolean; Run feature generation component')
def user_input(test_sequences, structure, features):
    '''
    Function for user to interact with.
    Test sequences is a two-column csv

    Args:
        TODO


    Returns:
        CSV with evaluation results.
    '''

    # user blast component
    """
    Input: csv file (list of sequences) + pairing and PDBids and UniprotIDs are optional
    Output:chunked csv for HMMR and structure components.
        Note: csv with specific parameters will be generated for specific component.
    Params:
    """
    logger = pairpro.utils.start_logger_if_necessary(
        LOGNAME, LOGFILE, LOGLEVEL, filemode='w')
    # convert csv to pandas dataframe
    df = pd.read_csv(test_sequences)

    # blast df has sequences and alignment metrics, PID that is unique for
    # each row
    df, con = make_blast_df(df)

    if len(df) < 1000:
        logger.info(
            'Running HMMER via the API as there are less than a 1000 sequences.')
        # subject_search = pairpro.hmmer.hmmerscanner(df, 'subject', 20, 20, HMMER_OUT_DIR)
        # query_search = pairpro.hmmer.hmmerscanner(df, 'query', 20, 20, HMMER_OUT_DIR)
        subject_scan = pairpro.hmmer.run_hmmerscanner(
            df, 'subject', 20, 20, HMMER_OUT_DIR)
        query_scan = pairpro.hmmer.run_hmmerscanner(
            df, 'query', 20, 20, HMMER_OUT_DIR)
        jaccard_threshold = 0.5
        # Get file pairs and calculate similarity for each pair
        file_pairs = pairpro.hmmer.get_file_pairs_API(HMMER_OUT_DIR)
        logger.info(
            f"Processing {len(file_pairs)} file pairs in {HMMER_OUT_DIR}")
        results = {}
        for file1, file2 in file_pairs:
            logger.info(f"Processing {file1} and {file2}")
            output_file = f"{PARSED_HMMER_OUT_DIR}functional_API_output.csv"
            similarity_scores = pairpro.hmmer.calculate_similarity_API(
                file1, file2, jaccard_threshold)
            pairpro.hmmer.write_function_output_API(
                similarity_scores, output_file)
            results[(file1, file2)] = output_file
        logger.info('Finished running HMMER via the API.')
    else:
        logger.info(
            'Running HMMER locally as there are more than a 1000 sequences.')
        logger.debug(f"Running HMMER locally with {len(df)} sequences.")
        chunk_size = 5000
        njobs = 4
        protein_chunks = [df[i:i + chunk_size]
                          for i in range(0, len(df), chunk_size)]
        logger.info(f'Running HMMER locally with {njobs} CPUs.')

        Parallel(
            n_jobs=njobs)(
            delayed(pairpro.hmmer.user_local_hmmer_wrapper_query)(
                chunk_index,
                PRESS_PATH,
                protein_chunks,
                HMMER_OUT_DIR) for chunk_index,
            protein_chunks in enumerate(protein_chunks))

        Parallel(
            n_jobs=njobs)(
            delayed(pairpro.hmmer.user_local_hmmer_wrapper_subject)(
                chunk_index,
                PRESS_PATH,
                protein_chunks,
                HMMER_OUT_DIR) for chunk_index,
            protein_chunks in enumerate(protein_chunks))

        file_pairs = pairpro.hmmer.get_file_pairs_user(HMMER_OUT_DIR)
        logger.info(
            f"Processing {len(file_pairs)} file pairs in {HMMER_OUT_DIR}")
        jaccard_threshold = 0.5
        results = {}
        for file1, file2 in file_pairs:
            logger.info(f"Processing {file1} and {file2}")
            file_chunk_index = int(file1.split("_")[-1].split(".")[0])
            output_file = f"{PARSED_HMMER_OUT_DIR}{file_chunk_index}_functional_output.csv"
            similarity_scores = pairpro.hmmer.calculate_similarity_user(
                file1, file2, jaccard_threshold)
            pairpro.hmmer.write_function_output_API(
                similarity_scores, output_file)
            results[(file1, file2)] = output_file
        logger.info('Finished running HMMER locally.')

        # checking if the parsed output is appended to table
        con.execute(
            """CREATE TABLE hmmer_results AS SELECT * FROM read_csv_auto('./data/user/parsed_hmmer_out/*.csv', HEADER=TRUE)""")
        con.execute(
            f"""ALTER TABLE protein_pairs ADD COLUMN hmmer_match BOOLEAN""")
        con.execute(f"""UPDATE protein_pairs AS f
        SET hmmer_match = hmmer.functional::BOOLEAN
        FROM hmmer_results AS hmmer
        WHERE hmmer.file1 = f.pair_id
        """)
        logger.info('Finished appending parsed HMMER output to table.')

    df = con.execute("""SELECT query, subject, bit_score, local_gap_compressed_percent_id,
    scaled_local_query_percent_id, scaled_local_symmetric_percent_id,
    query_align_len, query_align_cov, subject_align_len, subject_align_cov,
    LENGTH(query) AS query_len, LENGTH(subject) AS subject_len, hmmer_match FROM protein_pairs""").df()

    # ML component
    if structure:
        target = ['hmmer_match', 'structure_match']
    else:
        target = 'hmmer_match'

    # load model
    model = joblib.load(MODEL_PATH)

    # run model
    evaluate_model_wrapper(model, df, target, structure, features)

    logger.info('Finished running user input script.')


if __name__ == "__main__":
    logger = pairpro.utils.start_logger_if_necessary(
        LOGNAME, LOGFILE, LOGLEVEL, filemode='w')
    logger.info('Starting user input script.')

    # prepare output file
    try:
        os.makedirs(HMMER_OUT_DIR, exist_ok=True)
    except OSError as e:
        logger.error(f'Error creating directory: {e}')

    logger.info(f'Created HMMER output directory: {HMMER_OUT_DIR}')

    # prepare output file
    try:
        os.makedirs(PARSED_HMMER_OUT_DIR, exist_ok=True)
    except OSError as e:
        logger.error(f'Error creating directory: {e}')
    logger.info(
        f'Created parsed HMMER output directory: {PARSED_HMMER_OUT_DIR}')

    logger.info('Running user input script.')
    user_input()
