import sys
import logging
import os

# library dependencies
import pandas as pd
import joblib
from joblib import Parallel, delayed

# local dependencies
## machine learning
from PairPro.evaluate_model import evaluate_model
from PairPro.train_val_wrapper import train_val_wrapper
from PairPro.train_val_input_cleaning import columns_to_keep

#need to understand how to import the trained model from main
# from PairPro.main import train_model

## build DB
from PairPro.preprocessing import connect_db, build_fafsa
from PairPro.user_blast import make_blast_df

## hmmer
import PairPro.hmmer
import PairPro.utils


##structure
from PairPro.structure import download_structures, run_fatcat

### Paths
##ML Paths
MODEL_PATH = './data/models/'

def user_input(test_sequences, model, output_path:str):
    '''
    Function for user to interact with.
    Test sequences is a two-column csv

    Args:
        csv file (list of sequences)
        pairing and PDBids and UniprotIDs are optional
        

    Returns:
        CSV with evaluation results.
    '''

    # user blast component
    """
    Input: csv file (list of sequences) + pairing and PDBids and UniprotIDs are optional
    Output:chunked csv for HMMR and structure components.
        Note: csv with specific parameters will be generated for specific component.
    Params:
    """
    ## convert csv to pandas dataframe
    df = pd.read_csv(test_sequences)

    ## blast df has sequences and alignment metrics, PID that is unique for each row
    df = make_blast_df(df)

    # hmmer component
    """
    Input: Dataframe from user blast component
    Output: CSV (Meso_PID, Thermo_PID, Boolean)
    """
    user_boolean = make_target(df)

    ## developing method to generate single ID
    df = pd.merge(df, user_boolean, on=['ID'])

    # Structure component
    ## chau update code to return only boolean and append to the df
    download_structures(df, pdb_column, u_column, pdb_dir)
    df = run_fatcat(df, pdb_dir)

    #at this point, df is: user_blast + hmmer boolean + chau boolean

    # Machine Learning Component
    """
    Input: Dataframe that has been updated with user_blast + hmmer boolean + structure boolean
    Output: csv files with 6 columns (seq1, seq2, protein_match (Humood + Amin), protein_match (Chau) protein_match(ML))
    Params: Base environment + iFeatureOmega dependencies 
    """

    #make evaluation into four class classifier (neither true, hmmer true, structure true, both true)
    model = joblib.load(MODEL_PATH)
    evaluation = evaluate_model(df, model, output_path:str)
    
    return evaluation