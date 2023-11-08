import pandas as pd
import numpy as np
import torch
import esm
import os
import subprocess
import re

from pairpro.train_val_featuregen import get_fasta_from_dataframe
# from scripts.extract import main


def run_extract():
    # command for generating esm embeddings
    command = ['python', 'scripts/extract.py', 'esm2_t33_650M_UR50D', './notebooks/seq_a.fasta', 'data/pytorch_embeddings/', '--include', 'mean']

    # Run the command and capture the output
    output = subprocess.run(command, capture_output=True, text=True)

    # Print the captured output
    print("Output:", output.stdout)


## add embeddings to dataframe

def create_embeddings(OUTPUT_PATH='./data/pytorch_embeddings'):

    # Create an empty dataframe
    embeddings = pd.DataFrame(columns=['prot_pair_index', 'tensor'])

    # Iterate over the files in the directory
    for filename in os.listdir(OUTPUT_PATH):
        if filename.endswith('.pt'):  # Adjust the file extension as needed
            file_path = os.path.join(OUTPUT_PATH, filename)
            
            # Load the PyTorch file
            data = torch.load(file_path)
            
            # Extract the tensor(s) from the loaded data
            if isinstance(data, torch.Tensor):
                tensors = [data]
            elif isinstance(data, dict):
                tensors = list(data.values())
            else:
                # Handle other cases based on your specific file structure
                continue
            
            # Append the filename and tensor(s) to the dataframe
            embeddings = embeddings.append({'prot_pair_index': int(re.sub('.pt', '', filename)), 'tensor': tensors[1][33].numpy()}, ignore_index=True)

    #expand tensor values into individual columns
    embeddings = pd.concat([embeddings.drop('tensor', axis=1),
                    embeddings['tensor'].apply(lambda x: pd.Series(x))], axis=1)
    
    return embeddings
 
def update_dataframe(dataframe):

    """Generates embeddings from sequences in dataframe.

    Returns: Dataframe with ESM vector representation for each sequence.
        _type_: Pandas dataframe

    NOTE: Currently generates embeddings for meso seq only
          Could probably convert this to a script, need to figure
          out what fits into train_model better
    """

    #create fasta files
    get_fasta_from_dataframe(dataframe, 'seq_a.fasta', 'seq_b.fasta')

    #run extraction script
    run_extract()

    # create embeddings
    embeddings = create_embeddings()

    # merge two dataframes together
    dataframe = pd.merge(embeddings, dataframe, on='prot_pair_index')

    return dataframe