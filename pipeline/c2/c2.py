'''
This package samples ValidProt training data to represent taxa and protein pairs across the whole distribution of the ValidProt database.

Functions:
    validprot_sample: Generates a DataFrame containing a random or over/undersampled set from the ValidProt database.
    validprot_describe: Not supported in initial release. Provides analytics on ValidProt dataset to inform sampling parameters.
'''

# General strategy: create N-dimensional histogram where n is number of features worth sampling evenly. Function calculates the variance in height of these n-dimensional bins, then tries to minimize this variance iteratively. Potentially include a PCA 'fidelity' score to reduce dimensions and speed up calculations.
# How do we iterate? Take random sample of size = rate*size. Create n-d histogram. Calculate variance, mean, median, and min/max bin height. Reduce large bins to median height, increase small bins by duplicating samples, then cull to some tolerance % identical.
# Repeat sample with rate*size

import pandas as pd
import numpy as np

import time
import duckdb
import c0

def validprot_sample(df, stat_samp = 'random', exclude = None, plot_results = True, size, n_bins, tolerance, step_size):
    '''
    Generates sample DataFrame for downstream processing toward ValidProt training. Option to save analytics showing both sample
    and population distribution for each feature. When selecting 'random' sampling, the function simply passes df untouched, 
    assuming size has not changed. 
    
    Args:
        df (pandas.DataFrame): DataFrame formatted from component 1.
        stat_samp (str): One of ['random', 'boost_tails'].
        
    
    Returns:
    
    Raises:
    '''
    if stat_samp == 'random':
        
        sample_df = df.sample(size)
        
   
    if stat_samp = 'boost_tails':
        
        sample_df = df.sample(size)
    
    sample_df = pd.DataFrame()
    
    return sample_df

def validprot_describe(df, exclude):
    '''
    '''