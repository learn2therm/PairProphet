"""
TODO later
"""
# system dependencies
import os
import sys
import logging

# library dependencies
import duckdb as ddb
import numpy as np
from timeit import default_timer as timer

# local dependencies
import pairpro.utils as pp_utils

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

    # create directory for sample data
    try:
        os.makedirs('./data/OMA/samples', exist_ok=True)
    except OSError as e:
        logger.error(f'Error creating directory: {e}')
    logger.info(f'Created directory to store OMA samples: ./data/OMA/samples')

    # connect to database
    con = ddb.connect(OMA_DIR, read_only=False)
    logger.info("Connected to database %s", OMA_DIR)


    for sample in [100, 10000, 1000000]:
        logger.info(f'Starting sample {sample}')
        start_time = timer()
        df = pp_utils.oma_sample(con=con, size=sample)
        df.to_csv(f'./data/OMA/samples/oma_sample_{sample}.csv', index=False)
        logger.info(f'Completed sample {sample}: {timer()-start_time} seconds')

    # close connection
    con.close()
    logger.info("Closed connection to database %s", OMA_DIR)

