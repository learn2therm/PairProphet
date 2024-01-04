"""
NOTE: This script takes a long time to run.
The following is a test script to develop a method for downloading OMA data and building a database programmatically.

"""
# system dependencies
import datetime
from ftplib import FTP
import gzip
import logging
import os
import shutil

# library dependencies
import duckdb as ddb
from timeit import default_timer as timer

# local dependencies
import pairpro.utils as pp_utils

####################
### PATHS & VARS ###
####################
# paths
RAW_DATA_DIR = "./data/OMA"
DB_DIR = "./tmp"
DB_NAME = "oma.db"

# set up ftp
FTP_ADDRESS = 'omabrowser.org'

# variables
pair_file = '/All/oma-pairs.txt.gz'
seq_file = '/All/oma-seqs.fa.gz'
uniprot_file = '/All/oma-uniprot.txt.gz'
prok = '/All/prokaryotes.cdna.fa.gz'
# parquet dir?

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

def download_ftp_file(server, remote_file, local_file):
    """
    TODO: add docstring
    """
    ftp = FTP(server)
    ftp.login(user="anonymous", passwd='')
    ftp.cwd('/')  # Set the current directory (if needed)

    with open(local_file, 'wb') as file:
        ftp.retrbinary('RETR ' + remote_file, file.write)

    ftp.quit()

################
# Main script #
################



if __name__ == "__main__":

    # Initialize logger
    logger = pp_utils.start_logger_if_necessary(
        LOGNAME, LOGFILE, LOGLEVEL, filemode='w')
    logger.info("Starting script. Logging to %s", LOGFILE)

    # create directory for raw data
    try:
        os.makedirs(RAW_DATA_DIR, exist_ok=True)
    except OSError as e:
        logger.error(f'Error creating directory: {e}')

    logger.info(f'Created directory to store raw OMA data: {RAW_DATA_DIR}')

    # start timer
    time_start = timer()

    # download the raw data
    logger.info("Downloading raw data from OMA FTP server")

    logger.info("Downloading OMA uniprot file. Saving to %s", f"{RAW_DATA_DIR}/{os.path.basename(uniprot_file)}. Should take a few minutes.")
    download_ftp_file(FTP_ADDRESS, uniprot_file, f"{RAW_DATA_DIR}/{os.path.basename(uniprot_file)}")

    logger.info("Downloading OMA sequence file. Saving to %s", f"{RAW_DATA_DIR}/{os.path.basename(seq_file)}. Should take a couple of hours.")
    # continue here later


    # unzip the downloaded files
    logger.info("Unzipping downloaded files")

    # extracted_files? (list of extracted files?)

    with gzip.open(RAW_DATA_DIR, 'rb') as f_in:
        with open(extracted_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    logger.info(f'Extracted the file can be found in: {extracted_file}')

    # save metrics
    date_pulled = str(datetime.datetime.now().strftime("%m/%d/%Y"))
    with open('./data/pfam/pfam_pulled_timestamp', 'w') as file:
        file.write(date_pulled)
    logger.info(f'Wrote timestamp to file: {date_pulled}')

    