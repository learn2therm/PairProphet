"""
A quick test script to run the blast function and check the output.

Note: 
- This script is not intended to be run in production. It is only for testing purposes.
- The script can only run once you've generated the nuprok pairs from the OMA database 
    using the script `dev-OMA_nuprok.py`.
"""
# system dependencies
import os
import sys

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

# CLI arguments for sbatch
num_cpus = int(sys.argv[1])


if __name__ == "__main__":
    # initialize logger
    logger = pp_utils.start_logger_if_necessary(
        LOGNAME, LOGFILE, LOGLEVEL, filemode='w')
    logger.info("Starting script. Logging to %s", LOGFILE)

    # connect to database
    con = ddb.connect(OMA_DIR, read_only=False)
    logger.info("Connected to database %s", OMA_DIR)

    start_time = timer()
    logger.info(f"Obtaining combined nuprok pairs")
    df = con.execute("SELECT * FROM combined_pairs").fetchdf()

    logger.info(f"Starting blast function")
    blast_df, con = pp_up.make_blast_df(df, cpus=num_cpus, path=OMA_DIR)
    logger.info(f"Completed blast function: {timer()-start_time} seconds")

    # close connection
    con.close()
    logger.info("Closed connection to database %s", OMA_DIR)
    logger.info("Completed database operations!!")

    # save aligned dataframe
    logger.info("Saving aligned dataframe...")
    blast_df[['protein1_alphafold_id', 'protein2_alphafold_id']] = df[['protein1_uniprot_id', 'protein2_uniprot_id']]
    blast_df.to_csv('./data/OMA/samples/oma_sample_2m_aligned.csv', index=False)
    logger.info("Completed alignment. Saved to ./data/OMA/samples/oma_sample_2m_aligned.csv")