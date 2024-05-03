"""Ingest raw SIFTS data

Source
------

"""
# system dependecies
import datetime
import gzip
import logging
import os
import shutil
import tarfile


# library dependecies


# local dependencies
import pairpro.utils as pp_utils
import pairpro.structures as pp_struct

# get environmental variables
if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
    LOGLEVEL = getattr(logging, LOGLEVEL)
else:
    LOGLEVEL = logging.INFO
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'


# URL for a specific SIFTS file (e.g., UniProt to PDB mappings)
file_url = 'https://ftp.ebi.ac.uk/pub/databases/msd/sifts/csv/pdb_chain_uniprot.csv'
save_path = './data/SIFTS/pdb_chain_uniprot.csv'


if __name__ == "__main__":
    # set up logging
    logger = pp_utils.start_logger_if_necessary(
    LOGNAME, LOGFILE, LOGLEVEL, filemode='w')

    # get the logger in subprocesses
    # structure_logger = logging.getLogger('pairpro.structures')
    # structure_logger.setLevel(getattr(logging, LOGLEVEL))
    # structure_logger.addHandler(logging.FileHandler(LOGFILE))
    
    try:
        os.makedirs('./data/SIFTS', exist_ok=True)
    except OSError as e:
        logger.error(f'Error creating directory: {e}')

    logger.info(f'Created directory: ./data/SIFTS')

    logger.info(f'Downloading SIFTS file: {file_url}...')
    # download the raw data
    pp_struct.ProteinDownloader.download_sifts_file(file_url, save_path)
    logger.info("downloaded SIFTS file")


    # save metrics
    date_pulled = str(datetime.datetime.now().strftime("%m/%d/%Y"))
    with open('./data/SIFTS/SIFTS_pulled_timestamp', 'w') as file:
        file.write(date_pulled)
    logger.info(f'Wrote timestamp to file: {date_pulled}')