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
import threading

# library dependencies
import duckdb as ddb
import urllib.request
from urllib.error import URLError, HTTPError
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

# set up url address+files
# pair_file = '/All/oma-pairs.txt.gz
# seq_file = '/All/oma-seqs.fa.gz'
# uniprot_file = '/All/oma-uniprot.txt.gz'
# prok = '/All/prokaryotes.cdna.fa.gz'
File_urls = [
    'https://omabrowser.org/All/oma-uniprot.txt.gz', # uniprot mapping
    'https://omabrowser.org/All/oma-pairs.txt.gz', # pairwise orthologs
    'https://omabrowser.org/All/oma-seqs.fa.gz', # protein sequences
    'https://omabrowser.org/All/prokaryotes.cdna.fa.gz' # cDNA Prokaryotes
]
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

class FileDownloader:
    def __init__(self, logger):
        self.logger = logger
        self.download_directory = RAW_DATA_DIR
         # create directory for raw data
        try:
            os.makedirs(RAW_DATA_DIR, exist_ok=True)
        except OSError as e:
            self.logger.error(f'Error creating directory: {e}')
        self.logger.info(f'Created directory to store raw OMA data: {RAW_DATA_DIR}')

    def download_file(self, url):
        """
        TODO: add docstring
        """
        self.logger.info(f"Downloading {url}")
        file_name = url.split('/')[-1]
        destination = os.path.join(self.download_directory, file_name)
        try:
            with urllib.request.urlopen(url) as response, open(destination, 'wb') as out_file:
                shutil.copyfileobj(response, out_file)
                self.logger.info(f"Downloaded {url} to {destination}")
        except Exception as e:
            self.logger.error(f"Error downloading {url}: {e}")
            return None
        return destination

    def unzip_file(self, gz_path):
        """
        Unzip a .gz file to a destination path.
        """
        if gz_path:
            self.logger.info(f"Unzip file method called for {gz_path}")  # Test log
            dest_path = gz_path.rsplit('.gz', 1)[0]
            try:
                with gzip.open(gz_path, 'rb') as f_in:
                    with open(dest_path, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                self.logger.info(f"Unzipped {gz_path} to {dest_path}")
                return dest_path
            except Exception as e:
                self.logger.error(f"Error unzipping {gz_path}: {e}")
                return None
            
    def delete_file(self, file_path):
        """
        Delete a file from system.
        """
        if (file_path and os.path.exists(file_path)) and (file_path.endswith('.gz')):
            self.logger.info(f"Attempting to delete file: {file_path}")  # Log the path of the file being deleted
            try:
                os.remove(file_path)
                self.logger.info(f"Deleted {file_path}")
            except Exception as e:
                self.logger.error(f"Error deleting {file_path}: {e}")

    def omabrowser_download(self, url):
        """
        wrapper function for downloading files from OMA browser
        """
        self.logger.info(f"Starting download process for {url}")
        download_file = self.download_file(url)
        self.logger.info(f"Downloaded file path: {download_file}")

        if download_file:
            unzipped_file_path = self.unzip_file(download_file)
            self.logger.info(f"Unzipped file path: {unzipped_file_path}")

            self.delete_file(download_file)
            self.logger.info(f"Deleted file path: {download_file}")
        else:
            self.logger.error(f"Download failed for {url}, skipping unzip and delete.")
    
################
# Main script #
################



if __name__ == "__main__":

    # Initialize logger
    logger = pp_utils.start_logger_if_necessary(
        LOGNAME, LOGFILE, LOGLEVEL, filemode='w')
    logger.info("Starting script. Logging to %s", LOGFILE)

    # # Pass logger to FileDownloader class. Initialize FileDownloader
    # downloader = FileDownloader(logger)


    # threads = []
    # logger.info("Downloading raw data from OMA server")

    # time_start = timer()
    # for url in File_urls:
    #     thread = threading.Thread(target=downloader.omabrowser_download, args=(url,))
    #     thread.start()
    #     threads.append(thread)

    # for thread in threads:
    #     thread.join()

    # logger.info("Finished downloading raw data from OMA server. Time elapsed: %s", timer() - time_start)
    

    # # save metrics
    # date_pulled = str(datetime.datetime.now().strftime("%m/%d/%Y"))
    # with open('./data/OMA/oma_pulled_timestamp', 'w') as file:
    #     file.write(date_pulled)
    # logger.info(f'Wrote timestamp to file: {date_pulled}')

    try:
        os.makedirs(f"{RAW_DATA_DIR}/parquet", exist_ok=True)
    except OSError as e:
        logger.error(f"Error creating directory: {e}")
    logger.info(f"Created directory to store parquet files: {RAW_DATA_DIR}/parquet")

    # create database in tmp dir
    time_start = timer()
    con = ddb.connect(f"{DB_DIR}/{DB_NAME}", config={"allow_unsigned_extensions": True})
    logger.info(f"Created database at {DB_DIR}/{DB_NAME}")
    # install fasql extension
    logger.info("Installing fasql extension...")
    con.execute("SET custom_extension_repository='dbe.wheretrue.com/fasql/latest';")
    con.execute("INSTALL fasql;")
    con.execute("LOAD fasql;")

    # create parquet files
    logger.info("Creating parquet files...")
    # assigning column names
    pairs_cols = ["protein1", "protein2", "mult", "group"]
    uniprot_cols = ["oma_id", "uniprot_id"]
    # from utils, input file path as well as output parquet directory
    pp_utils.split_txt(in_path=f"{RAW_DATA_DIR}/oma-pairs.txt", out_dir=f"{RAW_DATA_DIR}/parquet", cols=pairs_cols)
    pp_utils.split_txt(in_path=f"{RAW_DATA_DIR}/oma-uniprot.txt", out_dir=f"{RAW_DATA_DIR}/parquet", cols=uniprot_cols)
    logger.info("Large text files converted to Parquet")

    # make pair SQL tables from parquet
    pp_utils.parqs_to_db(parq_dir=f"{RAW_DATA_DIR}/parquet/pairs", table_name=f"pairs", con=con, cols=['protein1', 'protein2'])
    pp_utils.parqs_to_db(parq_dir=f"{RAW_DATA_DIR}/parquet/pairs", table_name='uniprot', con=con)
    logger.info("SQL tables created within database. Parquet built within: %s", timer() - time_start)

    time_start2 = timer()
    # create/read fasta files
    con.execute(f"""CREATE OR REPLACE TABLE proteins AS SELECT description, sequence FROM read_fasta('{RAW_DATA_DIR}/oma-seqs.fa')""")
    logger.info("Built proteins in: %s", timer() - time_start2)

    time_start3 = timer()
    # create table with only prokrayotic pairs
    con.execute(f"""CREATE OR REPLACE TABLE prok_pairs AS SELECT
                protein1 AS protein1_oma_id,
                protein2 AS protein2_oma_id,
                FROM pairs 
                WHERE protein1 IN (SELECT description FROM read_fasta('{RAW_DATA_DIR}/prokaryotes.cdna.fa'))
                AND protein2 IN (SELECT description FROM read_fasta('{RAW_DATA_DIR}/prokaryotes.cdna.fa'))""")

    # Updates the prok_pairs table with uniprot ids. This and the previous steps must be split to avoid memory issues.
    con.execute("""CREATE OR REPLACE TABLE prok_pairs AS SELECT
                protein1_oma_id,
                protein2_oma_id,
                u1.uniprot_id AS protein1_uniprot_id,
                u2.uniprot_id AS protein2_uniprot_id,
                FROM prok_pairs
                LEFT JOIN uniprot u1 ON u1.oma_id = protein1_oma_id
                LEFT JOIN uniprot u2 ON u2.oma_id = protein2_oma_id""")
    
    logger.info("completed building OMA in: %s", timer() -time_start3)
    con.close()
    logger.debug("connection close. Script done.")

    