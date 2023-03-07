'''
This package samples ValidProt training data to represent taxa and protein pairs across the whole distribution of the ValidProt database.

Functions:
    validprot_sample: Generates a DataFrame containing a random or over/undersampled set from the ValidProt database.
'''

import pandas as pd
import numpy as np

import time
import duckdb
import c0

def validprot_sample(df, stat_samp = None, size):
    
    if stat_samp == None:
        
        sample_df = df.sample(size)
        return sample_df
   
    if stat_samp = representative:
        
        sample_df = df.sample(size)
    
    sample_df = pd.DataFrame()
    
    return sample_df

def explain_set(df, exclude):
    return