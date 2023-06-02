from Bio.PDB import PDBList
import os
import requests
import pandas as pd
import subprocess
import time
import tempfile

import asyncio
import httpx
import nest_asyncio

import duckdb as db
import numpy as np
import concurrent.futures
from joblib import Parallel, delayed

async def download_aff(session, url, filename):
    try:
        response = await session.get(url)
        if response.status_code == 200:
            with open(filename, 'wb') as f:
                f.write(response.content)
            print(f"Downloaded file: {filename}")
            return True
        else:
            print(f"Failed to download file: {filename}. Status code: {response.status_code}")
            return False
    except httpx.RequestError as e:
        print(f"Error while downloading file: {filename}. Exception: {str(e)}")
        return False

async def download_af(row, u_column, pdb_dir):
    uniprot_id = getattr(row, u_column)
    url = f'https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb'
    filename = f'{pdb_dir}/{uniprot_id}.pdb'

    async with httpx.AsyncClient() as client:
        success = await download_aff(client, url, filename)
        return success

def run_download_af_all(df, pdb_column, u_column, pdb_dir):
    nest_asyncio.apply()

    async def download_af_all():
        tasks = []
        success_count = 0

        if not os.path.exists(pdb_dir):
            os.makedirs(pdb_dir)

        for row in df.itertuples(index=False):
            if pd.isna(getattr(row, pdb_column)):
                task = asyncio.create_task(download_af(row, u_column, pdb_dir))
                tasks.append(task)

        results = await asyncio.gather(*tasks)
        success_count = sum(results)

        print(f"Successfully downloaded {success_count} files out of {len(df)}")

    asyncio.run(download_af_all())

def download_pdb(df, pdb_column, pdb_dir):
    pdbl = PDBList()
    for i, row in df.iterrows():
        pdb_id = row[pdb_column]
        if not pd.isna(pdb_id):  # check for NaN value in PDB IDs column
            pdbl.retrieve_pdb_file(pdb_id, pdir=pdb_dir, file_format='pdb')
            file_path = os.path.join(pdb_dir, f'pdb{pdb_id.lower()}.ent')
            if os.path.exists(file_path):
                os.rename(file_path, os.path.join(pdb_dir, f'{pdb_id}.pdb'))
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
    print(f"Execution time: {execution_time} seconds")
    pass

def compare_fatcat(p1_file, p2_file, pdb_dir, pair_id):
    cmd = ['FATCAT', '-p1', p1_file, '-p2', p2_file, '-i', pdb_dir, '-q']
    result = subprocess.run(cmd, capture_output=True, text=True)
    output = result.stdout
    p_value_line = next(line for line in output.split('\n') if line.startswith("P-value"))
    p_value = float(p_value_line.split()[1])
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
    p1_file = f'{p1}.pdb'
    p2_file = f'{p2}.pdb'
    if not os.path.exists(os.path.join(pdb_dir, p1_file)) or not os.path.exists(os.path.join(pdb_dir, p2_file)):
        return {'pair_id': row['pair_id'], 'p_value': np.nan}
    return compare_fatcat(p1_file, p2_file, pdb_dir, row['pair_id'])

def run_fatcat_dict_job(df, pdb_dir):
    start_time = time.time()
    p_values = []
    p_values = Parallel(n_jobs=-1)(delayed(process_row)(row, pdb_dir) for _, row in df.iterrows())
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Execution time: {execution_time} seconds")
    return p_values

if __name__ == '__main__':
    dff = pd.read_csv('../data/chau_test.csv')
    df = dff.sample(3)
    # Download the structures
    download_structure(df, 'meso_pdb', 'meso_pid', 'af')
    download_structure(df, 'thermo_pdb', 'thermo_pid', 'af')

    # Compare the structures and extract the p-values
    p_values = run_fatcat_dict_job(df, 'af')