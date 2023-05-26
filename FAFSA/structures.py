from Bio.PDB import PDBList
import os
import requests
import pandas as pd
import subprocess
import time
import tempfile

def download_structures(df, pdb_column, u_column, pdb_dir):
    start_time = time.time()  # Start measuring time
    pdbl = PDBList()
    if not os.path.exists(pdb_dir):
        os.makedirs(pdb_dir)
        
    for i, row in df.iterrows():
        pdb_id = row[pdb_column]
        uniprot_id = row[u_column]
        if not pd.isna(pdb_id):  # check for NaN value in PDB IDs column
            pdbl.retrieve_pdb_file(pdb_id, pdir=pdb_dir, file_format='pdb')
            file_path = os.path.join(pdb_dir, f'pdb{pdb_id.lower()}.ent')
            if os.path.exists(file_path):
                os.rename(os.path.join(file_path), os.path.join(pdb_dir, f'{pdb_id}.pdb'))
            else:
                pass
        elif isinstance(uniprot_id, str):  # download structure using UniProt ID
            url = f'https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb'
            response = requests.get(url)
            if response.ok:
                filename = f'{pdb_dir}/{uniprot_id}.pdb'
                with open(filename, 'wb') as f:
                    f.write(response.content)
                print(f"Downloaded file for {uniprot_id} to {filename}")
            else:
                print(f"Failed to download file for {uniprot_id}: {response.status_code} - {response.reason}")
        else:
            print(f"No PDB ID or UniProt ID available for index {i}")
        end_time = time.time()  # Stop measuring time
    execution_time = end_time - start_time
    print(f"Execution time: {execution_time} seconds")
    pass

def run_fatcat(df, pdb_dir):
    p_values = []  # List to store the extracted p-values
    rows_to_drop = []  # List to store the indices of rows to be dropped

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
            # Append the index of the row to the list of rows to be dropped
            rows_to_drop.append(index)
            continue

        # Set the FATCAT command and its arguments
        cmd = ['FATCAT', '-p1', p1_file, '-p2', p2_file, '-i', pdb_dir, '-q']
        
        # Run the FATCAT command and capture the output
        result = subprocess.run(cmd, capture_output=True, text=True)
        output = result.stdout

        # Find the line containing the p-value
        p_value_line = next(line for line in output.split('\n') if line.startswith("P-value"))

        # Extract the p-value
        p_value = float(p_value_line.split()[1])
        
        # Append the p-value to the list
        p_values.append(p_value)

    # Drop the rows with missing structure files from the dataframe
    df = df.drop(rows_to_drop)
    
    df.loc[:, 'p_value'] = p_values  # Use .loc to set the 'p_value' column
    return df