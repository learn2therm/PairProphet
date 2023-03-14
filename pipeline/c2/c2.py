'''
This package samples ValidProt training data to represent taxa and protein pairs across the whole distribution of the ValidProt database.

Functions:
    validprot_sample: Generates a DataFrame containing a random or over/undersampled set from the ValidProt database.
    validprot_describe: Not supported in initial release. Provides analytics on ValidProt dataset to inform sampling parameters.
    g_prime_i: Calculates g_prime vaues used in validprot_sample
'''

from typing import List
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

def validprot_sample(df, size, stat_samp = 'random', exclude = None, tol = 1e-04, density = 0.01):
    '''
    Generates sample DataFrame for downstream processing toward ValidProt training. Option to save analytics showing both sample
    and population distribution for each feature. When selecting 'random' sampling, the function simply passes df untouched, 
    assuming size has not changed.

    Args:
        df (pandas.DataFrame): DataFrame formatted from component 1.
        size (int): Size of output sample. Recommend << df for Frank Wolfe sampling.
        stat_samp (str): One of ['random', 'boost_tails'].
        exclude (list): List of non-numeric or unnecessary columns to be excluded for Frank Wolfe sampling. Leave as None for random.
        tol (float): D_lamb convergence condition for Frank Wolfe.
        density (float): Sampling density for Frank Wolfe.
        
    Returns:
        lamb (np.array): Probability distribution for optimal sampling of df.
        D_lamb_list (list): Saved list of D_lamb values per step.
        dd_lamb (list): Saved list of D_lamb changes per step. 

    Raises:
    '''
    if stat_samp == 'random':

        sample_df = df.sample(size)

    if stat_samp == 'frankwolfe':
        
        lamb, D_lamb_list, dd_lamb = FrankWolfe(df.to_numpy(), tol = tol)
        fw_idx = np.random.choice(range(df.shape[0]), size, p = lamb.flatten())
        sample_df = df.iloc[fw_idx]
                                               
    return sample_df, fw_idx

def FrankWolfe(X: np.ndarray, tol: float):
    '''
    Calculates probability distribution uing Frank Wolfe method to achieve D-optimal sampling of dataset.
    
    Args:
        X (np.array): Feature array. Must be numeric.
        tol (float): Convergence condition for D_lamb step. Saved as dd_lamb.
    '''
    D_lamb_list = []
    
    d = X.shape[1]
    n = X.shape[0]
    
    pulls = np.random.choice(list(range(n)), size=(2*d))
    values, counts_ = np.unique(pulls, return_counts=True)
    counts = []
    
    for i in range(n):
        
        if i in values:
            counts.append(counts_[np.argwhere(values==i)][0][0])
        else:
            counts.append(0)
            
    lamb = (np.array(counts)/(2*d)).reshape(-1,1)
    
    # start time after startup
    t = 2*d
    # the raw XXT matrix
    xs = []
    
    for x in X:
        xs.append(np.matmul(x.reshape(-1,1), x.reshape(1,-1)))
    XXT = np.array(xs)
    
    A = np.sum(lamb.reshape(-1,1,1)*XXT, axis=0)
    D_lamb = -np.log(np.linalg.det(A))
    D_lamb_list.append(D_lamb)
        
    while len(D_lamb_list) == 1:
        
        g_prime = np.concatenate([g_prime_i(X, A, i) for i in range(n)])
        It = np.argmin(g_prime)

        indicator = np.zeros((n,1))
        indicator[It] = 1
        
        # at this lambda and t are both at their end of previous iteration state, eg t has not been updated
        lamb = (lamb*t + indicator)/(t+1)
        t += 1
        
        A = np.sum(lamb.reshape(-1,1,1)*XXT, axis=0)
        D_lamb = -np.log(np.linalg.det(A))
        D_lamb_list.append(D_lamb)
        
    if len(D_lamb_list) > 1:
        while tol < abs(D_lamb-D_lamb_list[-2]):

            g_prime = np.concatenate([g_prime_i(X, A, i) for i in range(n)])
            It = np.argmin(g_prime)

            indicator = np.zeros((n,1))
            indicator[It] = 1

            # at this lambda and t are both at their end of previous iteration state, eg t has not been updated
            lamb = (lamb*t + indicator)/(t+1)
            t += 1

            A = np.sum(lamb.reshape(-1,1,1)*XXT, axis=0)
            D_lamb = -np.log(np.linalg.det(A))
            D_lamb_list.append(D_lamb)
    
    dd_lamb = [abs(D_lamb_list[i] - D_lamb_list[i-1]) for i in range(len(D_lamb_list)) if i > 0]
    return lamb, D_lamb_list, dd_lamb
    
def g_prime_i(X, A, i):
    '''
    Calculates g_prime for FrankWolfe function.
    
    Args:
        X (np.array): 
        A (float): 
        i (int): Index value of X
    Returns:
        g_prime (np.array): Calculated g_prime value for give i.  
    '''
    g_prime = -np.matmul(np.matmul(X[i].reshape(1,-1), np.linalg.inv(A)),  X[i].reshape(-1,1))
    return g_prime
