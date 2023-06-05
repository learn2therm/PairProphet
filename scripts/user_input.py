"""
This script is for user to interact with the pretrained model.

TODO:
- [ ] Add click functionality
"""

import sys
import logging
import os

# library dependencies
import pandas as pd
import joblib
from joblib import delayed, Parallel

# local dependencies
## machine learning
from pairpro.evaluate_model import evaluate_model
from pairpro.train_val_wrapper import train_val_wrapper
from pairpro.train_val_input_cleaning import columns_to_keep

#need to understand how to import the trained model from main
# from pairpro.main import train_model

## build DB
from pairpro.user_blast import make_blast_df

## hmmer
import pairpro.hmmer
import pairpro.utils


#structure
from pairpro.structures import download_structure, run_fatcat

### Paths
##ML Paths
MODEL_PATH = './data/models/'

## HMMER Paths
PRESS_PATH = './data/pfam/pfam'
HMMER_OUT_DIR = './data/user/hmmer_out'
PARSED_HMMER_OUT_DIR = './data/user/parsed_hmmer_out'

test_sequences = './data/50k_paired_seq.csv'


## get environmental variables
if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
    LOGLEVEL = getattr(logging, LOGLEVEL)
else:
    LOGLEVEL = logging.INFO
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'



def user_input(test_sequences, output_path:str, model=0):
    '''
    Function for user to interact with.
    Test sequences is a two-column csv

    Args:
        csv file (list of sequences)
        pairing and PDBids and UniprotIDs are optional
        

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
    logger = pairpro.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')
    ## convert csv to pandas dataframe
    df = pd.read_csv(test_sequences)

    ## blast df has sequences and alignment metrics, PID that is unique for each row
    df, con = make_blast_df(df)

    if df< 1000:
        logger.info('Running HMMER via the API as there are less than a 1000 sequences.')
        pairpro.hmmer.hmmerscanner(df, 20, 20, HMMER_OUT_DIR)
        pairpro.hmmer.run_hmmerscanner()
    else:
        logger.info('Running HMMER locally as there are more than a 1000 sequences.')
        chunk_size = 1000
        njobs = 4
        protein_chunks = [df[i:i + chunk_size] for i in range(0, len(df), chunk_size)]
        logger.info(f'Running HMMER locally with {njobs} CPUs.')
        Parallel(n_jobs=njobs)(delayed(pairpro.hmmer.user_local_hmmer_wrapper(
            chunk_index,
            PRESS_PATH,
            protein_chunks,
            HMMER_OUT_DIR) for chunk_index, 
            protein_chunks in enumerate(protein_chunks)))
        logger.info('Finished running HMMER locally.')
        jaccard_threshold = 0.5
        vector_size = 2
        logger.info('Parsing HMMER output.')
        pairpro.hmmer.process_pair_user(con, vector_size, jaccard_threshold, PARSED_HMMER_OUT_DIR)

        # checking if the parsed output is appended to table
        con.execute("""CREATE TABLE hmmer_results AS SELECT * FROM read_csv_auto('/data/user/parsed_hmmer_out/*.csv', HEADER=TRUE)""")
        con.execute(f"""ALTER TABLE proteins_pairs ADD COLUMN hmmer_match BOOLEAN""")
        con.execute(f"""UPDATE proteins_pairs AS f
        SET hmmer_match = hmmer.functional::BOOLEAN
        FROM hmmer_results AS hmmer
        WHERE hmmer.pair_id = f.pair_id
        """)
        logger.info('Finished appending parsed HMMER output to table.')
        


#     # hmmer component
#     """
#     Input: Dataframe from user blast component
#     Output: CSV (Meso_PID, Thermo_PID, Boolean)
#     """
#     user_boolean = make_target(df)

#     ## developing method to generate single ID
#     df = pd.merge(df, user_boolean, on=['ID'])

#     # Structure component
#     ## chau update code to return only boolean and append to the df
#     download_structures(df, pdb_column, u_column, pdb_dir)
#     df = run_fatcat(df, pdb_dir)

#     #at this point, df is: user_blast + hmmer boolean + chau boolean

#     # Machine Learning Component
#     """
#     Input: Dataframe that has been updated with user_blast + hmmer boolean + structure boolean
#     Output: csv files with 6 columns (seq1, seq2, protein_match (Humood + Amin), protein_match (Chau) protein_match(ML))
#     Params: Base environment + iFeatureOmega dependencies 
#     """

#     #make evaluation into four class classifier (neither true, hmmer true, structure true, both true)
#     model = joblib.load(MODEL_PATH)
#     evaluation = evaluate_model(df, model, output_path:str)
    
#     return evaluation

if __name__ == "__main__":
    logger = pairpro.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')
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

    user_input(test_sequences, output_path = 'no')