"""
This module takes in a pandas dataframe containing Uniprot IDs and PDB IDs, download the pdb files and run FATCAT for structural alignment purpose.
Returns a Boolean for structure similarity. If no structures were found for proteins, that pair is dropped in the output file.
"""
from Bio.PDB import PDBList
import os
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
            with open(self.cache_file, 'w') as f:
                json.dump({}, f)

    def compare_fatcat(self, args):
        p1_file, p2_file, pair_id = args

        # Check cache for result
        with open(self.cache_file, 'r') as f:
            cache = json.load(f)
            if pair_id in cache:
                logger.debug(f"Using cached result for pair_id: {pair_id}")
                return cache[pair_id]

        # If not in cache, process the pair
        cmd = ['FATCAT', '-p1', p1_file, '-p2', p2_file, '-i', self.pdb_dir, '-q']
        result = subprocess.run(cmd, capture_output=True, text=True)
        output = result.stdout
        p_value_line = next(line for line in output.split('\n') if line.startswith("P-value"))
        p_value = float(p_value_line.split()[1])
        result = {'pair_id': pair_id, 'p_value': True if p_value < 0.05 else False}

        # Save result to cache
        with open(self.cache_file, 'w') as f:
            cache[pair_id] = result
            json.dump(cache, f)

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
        p_values = []
        with Pool() as pool:
            args_list = [self.process_row(row) for _, row in df.iterrows() if self.process_row(row) is not None]
            results = pool.map(self.compare_fatcat, args_list)
            p_values = [p_value for p_value in results if p_value is not None]

        with open(output_file, 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=['pair_id', 'p_value'])
            writer.writeheader()
            writer.writerows(p_values)
        return output_file

        

# async def download_aff(session, url, filename, semaphore):
#     try:
#         async with semaphore:
#             response = await session.get(url)
#             if response.status_code == 200:
#                 with open(filename, 'wb') as f:
#                     f.write(response.content)
#                 logger.debug(f"Downloaded file: {filename}")
#                 return True
#             else:
#                 logger.debug(f"Failed to download file: {filename}. Status code: {response.status_code}")
#                 return False
#     except httpx.RequestError as e:
#         logger.info(f"Error while downloading file: {filename}. Exception: {str(e)}")
#         return False

# async def download_af(row, u_column, pdb_dir, semaphore):
#     uniprot_id = getattr(row, u_column)
#     url = f'https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb'
#     filename = f'{pdb_dir}/{uniprot_id}.pdb'

#     async with httpx.AsyncClient(verify=False) as client:  # Disable SSL certificate verification
#         success = await download_aff(client, url, filename, semaphore)
#         return success

# def run_download_af_all(df, pdb_column, u_column, pdb_dir):
#     nest_asyncio.apply()

#     async def download_af_all():
#         semaphore = asyncio.Semaphore(500) # Specify the maximum number of concurrent requests
#         tasks = []
#         success_count = 0

#         if not os.path.exists(pdb_dir):
#             os.makedirs(pdb_dir)

#         for row in df.itertuples(index=False):
#             if pd.isna(getattr(row, pdb_column)):
#                 task = asyncio.create_task(download_af(row, u_column, pdb_dir, semaphore))
#                 tasks.append(task)

#         results = await asyncio.gather(*tasks)
#         success_count = sum(results)

#         logger.info(f"Successfully downloaded {success_count} files out of {len(df)}")

#     asyncio.run(download_af_all())

# def download_pdb(df, pdb_column, pdb_dir):
#     pdbl = PDBList()
#     pdbs = df[pdb_column].dropna().unique()
#     for p in pdbs:
#         pdbl.retrieve_pdb_file(p, pdir=pdb_dir, file_format='pdb')
#         file_path = os.path.join(pdb_dir, f'pdb{p.lower()}.ent')
#         if os.path.exists(file_path):
#             os.rename(file_path, os.path.join(pdb_dir, f'{p}.pdb'))
#         else:
#             pass

# def download_structure(df, pdb_column, u_column, pdb_dir):
#     start_time = time.time()  # Start measuring time
#     if not os.path.exists(pdb_dir):
#         os.makedirs(pdb_dir)
#     download_pdb(df, pdb_column, pdb_dir)
#     run_download_af_all(df, pdb_column, u_column, pdb_dir)
#     end_time = time.time()  # Stop measuring time
#     execution_time = end_time - start_time
#     logger.info(f"Execution time: {execution_time} seconds")

# def compare_fatcat(p1_file, p2_file, pdb_dir, pair_id):
#     # Set the FATCAT command and its arguments
#     cmd = ['FATCAT', '-p1', p1_file, '-p2', p2_file, '-i', pdb_dir, '-q']

#     # Run the FATCAT command and capture the output
#     result = subprocess.run(cmd, capture_output=True, text=True)
#     output = result.stdout

#     # Find the line containing the p-value
#     p_value_line = next(line for line in output.split('\n') if line.startswith("P-value"))

#     # Extract the p-value and convert it to a numeric value
#     p_value = float(p_value_line.split()[1])

#     # Check if p-value is less than 0.05 and assign True or False accordingly
#     if p_value < 0.05:
#         return {'pair_id': pair_id, 'p_value': True}
#     else:
#         return {'pair_id': pair_id, 'p_value': False}

# def process_row(row, pdb_dir):
#     if not pd.isna(row['meso_pdb']):
#         p1 = row['meso_pdb']
#     else:
#         p1 = row['meso_pid']

#     if not pd.isna(row['thermo_pdb']):
#         p2 = row['thermo_pdb']
#     else:
#         p2 = row['thermo_pid']

#     # Check if the structure files exist in the 'checking' folder
#     p1_file = f'{p1}.pdb'
#     p2_file = f'{p2}.pdb'
#     if not os.path.exists(os.path.join(pdb_dir, p1_file)) or not os.path.exists(os.path.join(pdb_dir, p2_file)):
#         # Assign NaN as the p-value instead of dropping the row
#         return None

#     return p1_file, p2_file, pdb_dir, row['pair_id']

# def run_fatcat_dict_job(df, pdb_dir, file):
#     p_values = []

#     with Pool() as pool:
#         args_list = [process_row(row, pdb_dir) for _, row in df.iterrows() if process_row(row, pdb_dir) is not None]
#         results = pool.map(compare_fatcat, args_list)

#         p_values = [p_value for p_value in results if p_value is not None]

#     with open(file, 'w') as csvfile:
#         writer = csv.DictWriter(csvfile, fieldnames=['pair_id', 'p_value'])
#         writer.writeheader()
#         writer.writerows(p_values)
#     return file

# Uncomment the following lines to run the code as a script
# if __name__ == '__main__':
#     dff = pd.read_csv('../data/chau_test.csv')
#     df = dff.sample(5)
#     # Download the structures
#     download_structure(df, 'meso_pdb', 'meso_pid', 'af')
#     download_structure(df, 'thermo_pdb', 'thermo_pid', 'af')

#     # Compare the structures and extract the p-values
#     p_values = run_fatcat_dict_job(df, 'af', 'structure_alig.csv')

# Uncomment the following lines to run the class methods
# if __name__ == '__main__':
#     # set up logger
#     logging.basicConfig(level=logging.INFO)
#     logger = logging.getLogger(__name__)

#     dff = pd.read_csv('./data/chau_test.csv')
#     df = dff.sample(5)
#     downloader = ProteinDownloader(pdb_dir=STRUCTURE_DIR)
#     downloader.download_structure(df, pdb_column="meso_pdb", u_column="meso_pid")

#     processor = FatcatProcessor(pdb_dir=STRUCTURE_DIR)
#     processor.run_fatcat_dict_job(df, output_file=f"{STRUCTURE_DIR}output.csv")
