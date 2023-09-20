"""
This module takes in a pandas dataframe containing Uniprot IDs and PDB IDs, download the pdb files and run FATCAT for structural alignment purpose.
Returns a Boolean for structure similarity. If no structures were found for proteins, that pair is dropped in the output file.
"""
from Bio.PDB import PDBList
import os
from functools import partial
import pandas as pd
import subprocess
import time
import logging
import json

import asyncio
import httpx
import nest_asyncio

import duckdb as db
import numpy as np
import csv
from multiprocessing import Pool

logger = logging.getLogger(__name__)

# Structure Paths
STRUCTURE_DIR = './data/structures/'
STRUCTURE_OUTPUT_DIR = './data/protein_pairs/structures/'


class ProteinDownloader:
    """
    Class for downloading protein structures from PDB and AlphaFold.
    """
    def __init__(self, pdb_dir, alphafold_url_format="https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v4.pdb", max_concurrent_requests=500):
        self.pdb_dir = pdb_dir
        self.alphafold_url_format = alphafold_url_format
        self.max_concurrent_requests = max_concurrent_requests
        if not os.path.exists(self.pdb_dir):
            os.makedirs(self.pdb_dir)

    async def _download_aff(self, session, url, filename, semaphore):
        # Check if the file already exists (cache check)
        if os.path.exists(filename):
            logger.debug(f"File {filename} already exists. Skipping download.")
            return True

        try:
            async with semaphore:
                response = await session.get(url)
                if response.status_code == 200:
                    with open(filename, 'wb') as f:
                        f.write(response.content)
                    logger.debug(f"Downloaded file: {filename}")
                    return True
                else:
                    logger.debug(f"Failed to download file: {filename}. Status code: {response.status_code}")
                    return False
        except httpx.RequestError as e:
            logger.info(f"Error while downloading file: {filename}. Exception: {str(e)}")
            return False
        
    async def _download_af(self, row, u_column, semaphore):
        uniprot_id = getattr(row, u_column)
        url = self.alphafold_url_format.format(uniprot_id)
        filename = f'{self.pdb_dir}/{uniprot_id}.pdb'

        async with httpx.AsyncClient(verify=False) as client: # Disable SSL certificate verification
            success = await self._download_aff(client, url, filename, semaphore)
            return success
        
    def run_download_af_all(self, df, pdb_column, u_column):
        nest_asyncio.apply()

        async def download_af_all():
            semaphore = asyncio.Semaphore(self.max_concurrent_requests)
            tasks = []
            success_count = 0

            for row in df.itertuples(index=False):
                if pd.isna(getattr(row, pdb_column)):
                    task = asyncio.create_task(self._download_af(row, u_column, semaphore))
                    tasks.append(task)

            results = await asyncio.gather(*tasks)
            success_count = sum(results)

            logger.info(f"Successfully downloaded {success_count} files out of {len(df)}")

        asyncio.run(download_af_all())

    def download_pdb(self, df, pdb_column):
        pdbl = PDBList()
        pdbs = df[pdb_column].dropna().unique()
        for p in pdbs:
            file_path = os.path.join(self.pdb_dir, f'{p}.pdb')
            # Check if the file already exists (cache check)
            if os.path.exists(file_path):
                logger.debug(f"File {file_path} already exists. Skipping download.")
                continue
            pdbl.retrieve_pdb_file(p, pdir=self.pdb_dir, file_format='pdb')
            renamed_path = os.path.join(self.pdb_dir, f'pdb{p.lower()}.ent')
            if os.path.exists(renamed_path):
                os.rename(renamed_path, file_path)
            else:
                pass
    
    def download_structure(self, df, pdb_column, u_column):
        start_time = time.time()  # Start measuring time
        self.download_pdb(df, pdb_column)
        self.run_download_af_all(df, pdb_column, u_column)
        end_time = time.time()  # Stop measuring time
        execution_time = end_time - start_time
        logger.info(f"Execution time: {execution_time} seconds")


class FatcatProcessor:
    def __init__(self, pdb_dir, cache_file="./tmp/fatcat_cache.json"):
        self.pdb_dir = pdb_dir
        self.cache_file = cache_file
        # Check if cache file exists, if not, create an empty one
        if not os.path.exists(self.cache_file):
            open(self.cache_file, 'a').close()  # Just create an empty file if it doesn't exist.



    ## caching stuff
    def append_to_cache(self, result):
        if not result:
            logger.warning(f'Trying to append an empty result to the cache: {result}')
        else:
            logger.debug(f'Appending result to cache: {result}')
        # Check if the result is valid
        if not result or "pair_id" not in result or "p_value" not in result:
            logger.error(f"Invalid result to be cached: {result}")
            return
        # append to cache file
        with open(self.cache_file, 'a') as f:
            f.write(json.dumps(result) + "\n")

    def read_cache(self):
        # Check if cache file exists, if not, create an empty one
        if not os.path.exists(self.cache_file):
            return {}
        cache = {}
        if os.path.exists(self.cache_file):
            with open(self.cache_file, 'r') as f:
                for line in f:
                    line = line.strip() # remove whitespace
                    if not line or line == '{}': # skip empty lines or lines with '{}
                        continue
                    try:
                        item = json.loads(line)
                        pair_id = item.get("pair_id")
                        if pair_id is not None:
                            cache[pair_id] = item
                        else:
                            logger.warning(f"Unexpected cache entry found: {item}")
                    except json.JSONDecodeError:
                        logger.error(f"Failed to decode cache entry: {line.strip()}")
        return cache
    
    def clear_cache(self):
        cache = self.read_cache()
        with open(self.cache_file, 'w') as f:
            for pair_id, item in cache.items():
                f.write(json.dumps(item) + "\n")
    ## end caching stuff


    def compare_fatcat(self, args, cache):
        p1_file, p2_file, pair_id = args

        # Check if the pair is already in cache
        if str(pair_id) in cache:  # ensure the pair_id is a string, as JSON keys are always strings
            cached_item = cache[str(pair_id)]
            if "pair_id" in cached_item:
                logger.debug(f"Using cached result for pair_id: {pair_id}")
                return cached_item
            else:
                logger.warning(f"Cache entry for pair_id {pair_id} does not contain 'pair_id' key. Processing again.")

        # If not in cache, process the pair
        cmd = ['FATCAT', '-p1', p1_file, '-p2', p2_file, '-i', self.pdb_dir, '-q']
        try:
            logger.debug(f"Running FATCAT for pair_id: {pair_id} with cmd: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=1500)
            output = result.stdout
            p_value_line = next(line for line in output.split('\n') if line.startswith("P-value"))
            p_value = float(p_value_line.split()[1])
            result = {'pair_id': pair_id, 'p_value': True if p_value < 0.05 else False}
            # Save result to cache
            self.append_to_cache(result)
        except Exception as e:
            logger.error(f'Error while running FATCAT for pair_id: {pair_id}. Exception: {str(e)}')
            result = None

        return result


    def process_row(self, row):
        if not pd.isna(row['meso_pdb']):
            p1 = row['meso_pdb']
        else:
            p1 = row['meso_pid']

        if not pd.isna(row['thermo_pdb']):
            p2 = row['thermo_pdb']
        else:
            p2 = row['thermo_pid']

        p1_file = f'{p1}.pdb'
        p2_file = f'{p2}.pdb'
        if not os.path.exists(os.path.join(self.pdb_dir, p1_file)) or not os.path.exists(os.path.join(self.pdb_dir, p2_file)):
            return None

        return p1_file, p2_file, row['pair_id']

    def run_fatcat_dict_job(self, df, output_file):
        # check inital cache state
        with open(self.cache_file, 'r') as f:
            initial_cache_content = f.read()
        logger.debug(f"Initial cache content: {initial_cache_content}")

        # load cache at the beginning of the job
        cache = self.read_cache()

        num_cores = os.cpu_count() # get the number of available cores
        p_values = []
        with Pool(processes=num_cores) as pool: # use the dynamic number of processes
            args_list = [self.process_row(row) for _, row in df.iterrows() if self.process_row(row) is not None]
            chunk_size = max(1, len(args_list) // (2 * num_cores))  # Adjust this as needed
            results = pool.map(partial(self.compare_fatcat, cache=cache), args_list, chunksize=chunk_size)
            p_values = [p_value for p_value in results if p_value is not None]

        with open(output_file, 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=['pair_id', 'p_value'])
            writer.writeheader()
            writer.writerows(p_values)
        return output_file


# # Uncomment the following lines to run the class methods
# if __name__ == '__main__':
#     # set up logger
#     logger = logging.getLogger(__name__)
#     logger.setLevel(logging.DEBUG)
#     fh = logging.FileHandler('./logs/dev-structures.log', mode='w')
#     fh.setFormatter(logging.Formatter('%(filename)s - %(asctime)s %(levelname)-8s %(message)s'))
#     # conditional check to avoid duplicate handlers
#     if len(logger.handlers) == 0:
#         logger.addHandler(fh)
#     else:
#         logger.handlers[-1] = fh
    



#     dff = pd.read_csv('./data/chau_test.csv')
#     df = dff.sample(5)
#     downloader = ProteinDownloader(pdb_dir=STRUCTURE_DIR)
#     downloader.download_structure(df, pdb_column="meso_pdb", u_column="meso_pid")
#     downloader.download_structure(df, pdb_column="thermo_pdb", u_column="thermo_pid")

#     processor = FatcatProcessor(pdb_dir=STRUCTURE_DIR)
#     processor.run_fatcat_dict_job(df, output_file=f"{STRUCTURE_DIR}test_output.csv")
