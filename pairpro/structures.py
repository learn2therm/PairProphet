"""
This module takes in a pandas dataframe containing Uniprot IDs and PDB IDs, download the pdb files and run FATCAT for structural alignment purpose.
Returns a Boolean for structure similarity. If no structures were found for proteins, that pair is dropped in the output file.
"""
from Bio.PDB import PDBList
import os
import pandas as pd
import subprocess
import time

import asyncio
import httpx
import nest_asyncio

import duckdb as db
import numpy as np
import csv
from joblib import Parallel, delayed

async def download_aff(session, url, filename):
    """
    Downloads a file asynchronously using an HTTP session.

    Args:
        session (httpx.AsyncClient): An HTTP session for making requests.
        url (str): The URL of the file to download.
        filename (str): The name of the file to save.

    Returns:
        bool: True if the file is successfully downloaded, False otherwise.
    """
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
    """
    Downloads AlphaFold files for a given row asynchronously.

    Args:
        row (pd.Series): The row containing the data for the download.
        u_column (str): The column name for the UniProt ID.
        pdb_dir (str): The directory to save the downloaded files.

    Returns:
        bool: True if the download is successful, False otherwise.
    """
    uniprot_id = getattr(row, u_column)
    url = f'https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb'
    filename = f'{pdb_dir}/{uniprot_id}.pdb'

    async with httpx.AsyncClient(verify=False) as client:  # Disable SSL certificate verification
        success = await download_aff(client, url, filename)
        return success

def run_download_af_all(df, pdb_column, u_column, pdb_dir):
    """
    Runs the asynchronous download of AlphaFold files for all rows in a DataFrame.

    Args:
        df (pd.DataFrame): The DataFrame containing the data for the downloads.
        pdb_column (str): The column name for the PDB ID.
        u_column (str): The column name for the UniProt ID.
        pdb_dir (str): The directory to save the downloaded files.
    
    Returns:
        files: pdb files containing structural information.
    """
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
    """
    Downloads PDB files for the given DataFrame based on PDB IDs.

    Args:
        df (pd.DataFrame): The DataFrame containing the PDB IDs.
        pdb_column (str): The column name for the PDB ID.
        pdb_dir (str): The directory to save the downloaded files.
    
    Returns: pdb files containing structural information.
    """
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
    """
    Downloads structure files for a DataFrame using AlphaFold and PDB.

    Args:
        df (pd.DataFrame): The DataFrame containing the data for the downloads.
        pdb_column (str): The column name for the PDB ID.
        u_column (str): The column name for the UniProt ID.
        pdb_dir (str): The directory to save the downloaded files.

    Returns: pdb files containing structural information.
    """
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
    """
    Compares two protein structures using FATCAT.

    Args:
        p1_file (str): The path to the first protein structure file.
        p2_file (str): The path to the second protein structure file.
        pdb_dir (str): The directory containing the structure files.
        pair_id (str): The ID of the protein pair.

    Returns:
        dict: A dictionary containing the pair ID and the p-value.
    """
    # Set the FATCAT command and its arguments
    cmd = ['FATCAT', '-p1', p1_file, '-p2', p2_file, '-i', pdb_dir, '-q']

    # Run the FATCAT command and capture the output
    result = subprocess.run(cmd, capture_output=True, text=True)
    output = result.stdout

    # Find the line containing the p-value
    p_value_line = next(line for line in output.split('\n') if line.startswith("P-value"))

    # Extract the p-value and convert it to a numeric value
    p_value = float(p_value_line.split()[1])

    # Check if p-value is less than 0.05 and assign True or False
    if p_value < 0.05:
        return {'pair_id': pair_id, 'p_value': True}
    else:
        return {'pair_id': pair_id, 'p_value': False}

def process_row(row, pdb_dir):
    """
    Processes a row of a DataFrame to compare protein structures using FATCAT.

    Args:
        row (pd.Series): The row containing the data for the comparison.
        pdb_dir (str): The directory containing the structure files.

    Returns:
        dict: A dictionary containing the pair ID and the p-value.
    """
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

def run_fatcat_dict_job(df, pdb_dir, file):
    """
    Runs the FATCAT comparison job on a DataFrame and saves the results to a file.

    Args:
        df (pd.DataFrame): The DataFrame containing the data for the comparison.
        pdb_dir (str): The directory containing the structure files.
        file (str): The path to the output file.

    Returns:
        csv: a csv file containing the pair ID and the p-value.
    """
    p_values = []  # List to store the extracted p-values

    # Parallelize the execution of the function using joblib
    p_values = Parallel(n_jobs=-1)(delayed(process_row)(row, pdb_dir) for _, row in df.iterrows())

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