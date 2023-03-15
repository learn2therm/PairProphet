'''
This package samples ValidProt training data to represent taxa and protein pairs across the whole distribution of the ValidProt database.

Functions:
    validprot_sample: Generates a DataFrame containing a random or over/undersampled set from the ValidProt database.
    validprot_describe: Not supported in initial release. Provides analytics on ValidProt dataset to inform sampling parameters.
    g_prime_i: Calculates g_prime vaues used in validprot_sample
'''

import numpy as np
import pandas as pd

def validprot_sample(df, size, stat_samp = 'random', tol = 1e-04):
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
        sample_df (pd.DataFrame): Sample of df generated using sampling function.
        fw_idx (list): Row indices selected by sampling function. Useful if df is derived from a DataFrame with identical rows and 
        more columns.

    Raises:
        AttributeError: Necessary features for ValidProt training must exist in input DataFrame.
    '''
    features = ['local_gap_compressed_percent_id', 'scaled_local_query_percent_id', 
        'scaled_local_symmetric_percent_id', 'query_align_len', 'subject_align_len', 
            'subject_align_cov', 'bit_score', 'm_protein_len', 't_protein_len']
    
    if False in [i in list(df.columns) for i in features]:
        raise AttributeError("""DataFrame must contain ValidProt features: local_gap_compressed_percent_id,
                             scaled_local_query_percent_id, scaled_local_symmetric_percent_id, query_align_len,
                             subject_align_len, subject_align_cov, bit_score, m_protein_len, t_protein_len]""")
        
    # Simply takes a random sample of input data and saves indices of the random choice.
    if stat_samp == 'random':
        
        fw_idx = pd.Index(np.random.choice(df.index, size))
        sample_df = df.iloc[fw_idx]

    # Applies Frank-Wolfe algorithm for D-optimal sampling. Recommend size << df rows.
    if stat_samp == 'frankwolfe':
        
        fw_array = df[features].to_numpy()
        lamb, D_lamb_list, dd_lamb = FrankWolfe(fw_array, tol = tol)
        fw_idx = np.random.choice(range(df.shape[0]), size, p = lamb.flatten())
        sample_df = df.iloc[fw_idx]
                                                  
    return sample_df, fw_idx


def FrankWolfe(X: np.ndarray, tol: float):
    '''
    Calculates probability distribution uing Frank Wolfe method to achieve D-optimal sampling of dataset.
    
    Args:
        X (np.array): Feature array. Must be numeric.
        tol (float): Convergence condition for D_lamb step. Saved as dd_lamb.
    
    Returns:
        lamb (np.array): Probability distribution solved by convergence of Frank-Wolfe calculation.
        D_lamb_list (list): List of D_lamb calculated at each iteration.
        dd_lamb (list): List of absolute delta D_lamb between each step. Final value should be below tol.
    '''
    
    D_lamb_list = []
    
    d = X.shape[1]
    n = X.shape[0]
    
    # Takes initial sample for first step of FW algorithm
    pulls = np.random.choice(list(range(n)), size=(2*d))
    values, counts_ = np.unique(pulls, return_counts=True)
    counts = []
    
    for i in range(n):
        
        if i in values:
            counts.append(counts_[np.argwhere(values==i)][0][0])
        else:
            counts.append(0)
    
    # Calculates initial probability distribution from initialization step
    lamb = (np.array(counts)/(2*d)).reshape(-1,1)
    
    t = 2*d
    xs = []
    
    # Calculate XXT valaues
    for x in X:
        xs.append(np.matmul(x.reshape(-1,1), x.reshape(1,-1)))
    XXT = np.array(xs)
    
    # Get more initial values
    A = np.sum(lamb.reshape(-1,1,1)*XXT, axis=0)
    D_lamb = -np.log(np.linalg.det(A))
    D_lamb_list.append(D_lamb)
        
    # Manages first step of iterator
    while len(D_lamb_list) == 1:
        
        g_prime = np.concatenate([g_prime_i(X, A, i) for i in range(n)])
        It = np.argmin(g_prime)

        indicator = np.zeros((n,1))
        indicator[It] = 1
        
        lamb = (lamb*t + indicator)/(t+1)
        t += 1
        
        A = np.sum(lamb.reshape(-1,1,1)*XXT, axis=0)
        D_lamb = -np.log(np.linalg.det(A))
        D_lamb_list.append(D_lamb)
        
    # Iterates FW until change in D_lamb is less than tol
    if len(D_lamb_list) > 1:
        while tol < abs(D_lamb-D_lamb_list[-2]):

            g_prime = np.concatenate([g_prime_i(X, A, i) for i in range(n)])
            It = np.argmin(g_prime)

            indicator = np.zeros((n,1))
            indicator[It] = 1

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
