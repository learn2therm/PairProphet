"""
TODO
"""
# system dependencies
import os
import sys
import logging

# library dependencies
import duckdb as ddb
import pandas as pd
from timeit import default_timer as timer

# local dependencies
import pairpro.utils as pp_utils
import pairpro.user_blast as pp_up

####################
### PATHS & VARS ###
####################
# paths
OMA_DIR = "./tmp/oma.db"

# venv variables
if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
else:
    LOGLEVEL = 'DEBUG' # change to INFO for production
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'


if __name__ == "__main__":
    # initialize logger
    logger = pp_utils.start_logger_if_necessary(
        LOGNAME, LOGFILE, LOGLEVEL, filemode='w')
    logger.info("Starting script. Logging to %s", LOGFILE)

    # connect to database
    con = ddb.connect(OMA_DIR, read_only=False)
    logger.info("Connected to database %s", OMA_DIR)

    start_time = timer()
    logger.info(f'Starting sample to generate nuprok pairs; 1.5M pairs')
    df = pp_utils.nuprok_sample(con=con, size=2000000)
    logger.info(f'Completed nuprok pairs: {timer()-start_time} seconds')

    logger.info(f"original dataframe: {df.shape}. Proceeding to process and clean the table")

    # create 'clean' nuprok pairs
    start_time = timer()
    con.execute("""CREATE OR REPLACE TABLE clean_nuprok_pairs 
                AS SELECT * FROM nuprok_pairs 
                WHERE protein1_uniprot_id <> protein2_uniprot_id AND protein1_sequence <> protein2_sequence 
                AND protein1_uniprot_id IS NOT NULL 
                AND protein2_uniprot_id IS NOT NULL 
                AND protein1_sequence IS NOT NULL 
                AND protein2_sequence IS NOT NULL""")
    logger.info(f'Completed processing nuprok pairs: {timer()-start_time} seconds')
    logger.info(f"clean_nuprok_pairs: {con.table('clean_nuprok_pairs').shape()}")

    # create 'bad' sample from clean_nuprok_pairs
    logger.info(f"creating 'bad' sample from clean_nuprok_pairs")
    start_time = timer()
    con.execute("""SELECT 'bad_pair_' || s1.rn AS pair_id, s1.protein1_uniprot_id, s2.protein2_uniprot_id, s1.protein1_sequence, s2.protein2_sequence 
                FROM (
                SELECT row_number() OVER (ORDER BY random()) AS rn, protein1_uniprot_id, protein1_sequence 
                FROM clean_nuprok_pairs 
                USING SAMPLE 50%) s1 
                JOIN (
                SELECT row_number() OVER (ORDER BY random()) AS rn, protein2_uniprot_id, protein2_sequence 
                FROM clean_nuprok_pairs 
                USING SAMPLE 50%) s2 ON s1.rn = s2.rn;""").df()
    
    


    # logger.info("Beginning to processing nuprok pairs")
    # logger.info(f"original dataframe: {df.shape}")

    # # drop duplicates IDs
    # logger.info("Dropping duplicates ids in nuprok pairs")
    # id_id = df[df.protein1_uniprot_id == df.protein2_uniprot_id].index
    # df.drop(id_id, inplace=True)
    # df.reset_index(drop=True, inplace=True)

    # # drop duplicates sequences
    # logger.info("Dropping duplicates sequences in nuprok pairs")
    # seq_id = df[df.protein1_sequence == df.protein2_sequence].index
    # df.drop(seq_id, inplace=True)
    # df.reset_index(drop=True, inplace=True)

    # # drop na
    # logger.info("Dropping na in nuprok pairs")
    # df.dropna(inplace=True)
    # df.reset_index(drop=True, inplace=True)
    # logger.info(f"processed dataframe: {df.shape}")

    # # close connection
    # con.close()
    # logger.info("Closed connection to database %s", OMA_DIR)

    # module_path=os.getcwd()

    # pp_up.__path__=module_path
    # logger.debug(f"CHECK-> module path: {module_path}")

    # logger.info("Starting to align nuprok pairs")
    # start = timer()
    # df_align, con = pp_up.make_blast_df(df[['protein1_sequence', 'protein2_sequence']], cpus=4, path='./data/OMA/samples/oma_sample_1.5m.db', module_path=module_path)
    # time_delta = timer() - start
    # logger.info(f'Alignment time: {time_delta} seconds')

    # # close connection
    # con.close()
    # logger.info("Closed connection to database %s", './data/OMA/samples/oma_sample_1.5m.db')

    # # save aligned dataframe
    # logger.info("Saving aligned dataframe...")
    # df_align[['protein1_alphafold_id', 'protein2_alphafold_id']] = df[['protein1_uniprot_id', 'protein2_uniprot_id']]
    # df_align.to_csv('./data/OMA/samples/oma_sample_1.5m_aligned.csv', index=False)
    # logger.info("Completed alignment. Saved to ./data/OMA/samples/oma_sample_1.5m_aligned.csv")