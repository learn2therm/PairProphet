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

import asyncio
import httpx
import nest_asyncio

import duckdb as db
import numpy as np
import csv
from joblib import Parallel, delayed

logger = logging.getLogger(__name__)

async def download_aff(session, url, filename, semaphore):
    try:
        async with semaphore:
            response = await session.get(url)
            if response.status_code == 200:
                with open(filename, 'wb') as f:
                    f.write(response.content)
                logger.info(f"Downloaded file: {filename}")
                return True
            else:
                logger.info(f"Failed to download file: {filename}. Status code: {response.status_code}")
                return False
    except httpx.RequestError as e:
        logger.info(f"Error while downloading file: {filename}. Exception: {str(e)}")
        return False

async def download_af(row, u_column, pdb_dir, semaphore):
    uniprot_id = getattr(row, u_column)
    url = f'https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb'
    filename = f'{pdb_dir}/{uniprot_id}.pdb'

    async with httpx.AsyncClient(verify=False) as client:  # Disable SSL certificate verification
        success = await download_aff(client, url, filename, semaphore)
        return success

def run_download_af_all(df, pdb_column, u_column, pdb_dir):
    nest_asyncio.apply()

    async def download_af_all():
        semaphore = asyncio.Semaphore(500) # Specify the maximum number of concurrent requests
        tasks = []
        success_count = 0

        if not os.path.exists(pdb_dir):
            os.makedirs(pdb_dir)

        for row in df.itertuples(index=False):
            if pd.isna(getattr(row, pdb_column)):
                task = asyncio.create_task(download_af(row, u_column, pdb_dir, semaphore))
                tasks.append(task)

        results = await asyncio.gather(*tasks)
        success_count = sum(results)

        logger.info(f"Successfully downloaded {success_count} files out of {len(df)}")

    asyncio.run(download_af_all())

def download_pdb(df, pdb_column, pdb_dir):
    pdbl = PDBList()
    pdbs = df[pdb_column].dropna().unique()
    for p in pdbs:
        pdbl.retrieve_pdb_file(p, pdir=pdb_dir, file_format='pdb')
        file_path = os.path.join(pdb_dir, f'pdb{p.lower()}.ent')
        if os.path.exists(file_path):
            os.rename(file_path, os.path.join(pdb_dir, f'{p}.pdb'))
        else:
            pass

def download_structure(df, pdb_column, u_column, pdb_dir):
    start_time = time.time()  # Start measuring time
    if not os.path.exists(pdb_dir):
        os.makedirs(pdb_dir)
    download_pdb(df, pdb_column, pdb_dir)
    run_download_af_all(df, pdb_column, u_column, pdb_dir)
    end_time = time.time()  # Stop measuring time
    execution_time = end_time - start_time
    logger.info(f"Execution time: {execution_time} seconds")

def compare_fatcat(p1_file, p2_file, pdb_dir, pair_id):
    # Set the FATCAT command and its arguments
    cmd = ['FATCAT', '-p1', p1_file, '-p2', p2_file, '-i', pdb_dir, '-q']

    # Run the FATCAT command and capture the output
    result = subprocess.run(cmd, capture_output=True, text=True)
    output = result.stdout

    # Find the line containing the p-value
    p_value_line = next(line for line in output.split('\n') if line.startswith("P-value"))

    # Extract the p-value and convert it to a numeric value
    p_value = float(p_value_line.split()[1])

    # Check if p-value is less than 0.05 and assign 1 or 0 accordingly
    if p_value < 0.05:
        return {'pair_id': pair_id, 'p_value': 1}
    else:
        return {'pair_id': pair_id, 'p_value': 0}

def process_row(row, pdb_dir):
    if not pd.isna(row['meso_pdb']):
        p1 = row['meso_pdb']
    else:
        p1 = row['meso_pid']

    if not pd.isna(row['thermo_pdb']):
        p2 = row['thermo_pdb']
    else:
        p2 = row['thermo_pid']

    # Check if the structure files exist in the 'checking' folder
    p1_file = f'{p1}.pdb'
    p2_file = f'{p2}.pdb'
    if not os.path.exists(os.path.join(pdb_dir, p1_file)) or not os.path.exists(os.path.join(pdb_dir, p2_file)):
        # Assign NaN as the p-value instead of dropping the row
        return None

    return compare_fatcat(p1_file, p2_file, pdb_dir, row['pair_id'])

def run_fatcat_dict_job(df, pdb_dir, njobs, file):
    p_values = []  # List to store the extracted p-values

    # Parallelize the execution of the function using joblib
    p_values = Parallel(n_jobs=njobs)(delayed(process_row)(row, pdb_dir) for _, row in df.iterrows())
    # p_values = Parallel(n_jobs=num_cores)(delayed(process_row)(row, pdb_dir) for _, row in df.iterrows())
    # Filter out None values
    p_values = [p_value for p_value in p_values if p_value is not None]

    with open(file, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=['pair_id', 'p_value'])
        writer.writeheader()
        writer.writerows(p_values)
    return file

# Uncomment the following lines to run the code as a script
# if __name__ == '__main__':
#     dff = pd.read_csv('../data/chau_test.csv')
#     df = dff.sample(5)
#     # Download the structures
#     download_structure(df, 'meso_pdb', 'meso_pid', 'af')
#     download_structure(df, 'thermo_pdb', 'thermo_pid', 'af')

#     # Compare the structures and extract the p-values
#     p_values = run_fatcat_dict_job(df, 'af', 'structure_alig.csv')