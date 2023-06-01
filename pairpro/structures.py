"""
Contains functions to download structures from PDB or AlphaFold2 and run FATCAT on them.
"""

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

async def download_aff(session, url, filename):
    """
    Downloads a file from the given URL using the provided session object and saves it with the specified filename.

    Args:
        session (httpx.AsyncClient): An asynchronous HTTP client session.
        url (str): The URL from which to download the file.
        filename (str): The name of the file to save the downloaded content.

    Returns:
        bool: True if the file was downloaded successfully, False otherwise.

    Raises:
        httpx.RequestError: If an error occurs during the HTTP request.
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
    Downloads an AlphaFold PDB file corresponding to the given row and saves it in the specified directory.

    Args:
        row (pd.Series): A row from a pandas DataFrame containing the necessary information.
        u_column (str): The name of the column in the DataFrame containing the UniProt ID.
        pdb_dir (str): The directory where the PDB file should be saved.

    Returns:
        bool: True if the file was downloaded successfully, False otherwise.
    """
    uniprot_id = getattr(row, u_column)
    url = f'https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb'
    filename = f'{pdb_dir}/{uniprot_id}.pdb'

    async with httpx.AsyncClient() as client:
        success = await download_aff(client, url, filename)
        return success

def run_download_af_all(df, u_column, pdb_dir):
    """
    Runs the asynchronous downloading of AlphaFold pdb files for all rows in the given DataFrame.

    Args:
        df (pd.DataFrame): The DataFrame containing the data for downloading PDB files.
        u_column (str): The name of the column in the DataFrame containing the UniProt IDs.
        pdb_dir (str): The directory where the PDB files should be saved.

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
    Downloads pdb files corresponding to the PDB IDs in the specified column of the DataFrame and saves them in the
    specified directory.

    Args:
        df (pd.DataFrame): The DataFrame containing the data for downloading PDB files.
        pdb_column (str): The name of the column in the DataFrame containing the PDB IDs.
        pdb_dir (str): The directory where the PDB files should be saved.

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
    Downloads PDB files and AlphaFold models for the given DataFrame and saves them in the specified directory. This
    function combines the functionality of `download_pdb` and `run_download_af_all`.

    Args:
        df (pd.DataFrame): The DataFrame containing the data for downloading PDB files and AlphaFold models.
        pdb_column (str): The name of the column in the DataFrame containing the PDB IDs.
        u_column (str): The name of the column in the DataFrame containing the UniProt IDs.
        pdb_dir (str): The directory where the PDB files and AlphaFold models should be saved.
    
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

def run_fatcat(df, pdb_dir):
    """
    Runs the FATCAT protein structure alignment tool on the PDB files in the specified directory and extracts the
    p-values. It updates the DataFrame with the calculated p-values.

    Args:
        df (pd.DataFrame): The DataFrame containing the data for running the FATCAT tool.
        pdb_dir (str): The directory where the PDB files are located.

    Returns:
        pd.DataFrame: The updated DataFrame with the added 'p_value' column.
        p_values: 0 False pair
                  1 True pair
                  NaN if structure file does not exist
    """
    p_values = []  # List to store the extracted p-values

    for index, row in df.iterrows():
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
            # Assign a p-value of 2 to the row instead of dropping it
            p_values.append(np.nan)
            continue

        # Set the FATCAT command and its arguments
        cmd = ['FATCAT', '-p1', p1_file, '-p2', p2_file, '-i', pdb_dir, '-q']
        
        # Run the FATCAT command and capture the output
        result = subprocess.run(cmd, capture_output=True, text=True)
        output = result.stdout

        # Find the line containing the p-value
        p_value_line = next(line for line in output.split('\n') if line.startswith("P-value"))

        # Extract the p-value and convert it to numeric value
        p_value = float(p_value_line.split()[1])
        
        # Check if p-value is less than 0.05 and assign 1 or 0 accordingly
        if p_value < 0.05:
            p_values.append(0)
        else:
            p_values.append(1)

    df.loc[:, 'p_value'] = p_values  # Use .loc to set the 'p_value' column
    return df

# if __name__ == "__main__":
#     # Read the dataframe
#     df = pd.read_csv('../data/pair_sample.csv')
#     df_sample = df.sample(5)
#     # Download structures from PDB or AlphaFold2
#     download_structure(df_sample, pdb_column='meso_pdb', u_column='meso_pid', pdb_dir='structures')
#     download_structure(df_sample, pdb_column='thermo_pdb', u_column='thermo_pid', pdb_dir='structures')
#     # Run FATCAT on the downloaded structures
#     df_result = run_fatcat(df_sample, pdb_dir='structures')