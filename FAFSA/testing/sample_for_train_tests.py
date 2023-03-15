'''
Unit testing script for FAFSA sample_for_train.
'''

import unittest
import os
import pandas as pd
import numpy as np

from FAFSA.sample_for_train import *

def get_clean_df(filename):
    
    features = ['local_gap_compressed_percent_id', 'scaled_local_query_percent_id', 
            'scaled_local_symmetric_percent_id', 'query_align_len', 'subject_align_len', 
                'subject_align_cov', 'bit_score', 'm_protein_len', 't_protein_len']
    
    # Get path for test dataset import
    db_path = os.path.abspath(os.path.join('..', 'data', filename))

    if os.path.exists(db_path) is False:
        raise ValueError(f'Could not find {filename} in current directory')
        
    df = pd.read_csv(db_path)
    
    return df[features]


class TestSample(unittest.TestCase):
    ''' 
    Tests for the fafsa_sample function.
    '''
    
    def test_smoke(self):
        '''
        Smoke test to ensure function runs both modes without error.
        '''
        assert c2.fafsa_sample(get_clean_df('learn2therm_sample_50k.zip').sample(100), 
                                   size = 10, stat_samp = 'frankwolfe')
        
    def test_oneshot_random(self):
        '''
        Test that correct sample data is returned from known inputs.
        '''
        test_df = get_clean_df('learn2therm_sample_50k.zip').sample(100).reset_index()
        df, idx = c2.fafsa_sample(test_df, size = 10, stat_samp = 'random')
        
        assert df.shape[0] == 10
        assert len(idx) == 10

        
    def test_oneshot_FW(self):
        '''
        Test that correct sample data is returned from known inputs.
        '''
        test_df = get_clean_df('learn2therm_sample_50k.zip').sample(100).reset_index()
        df, idx = c2.fafsa_sample(test_df, size = 10, stat_samp = 'frankwolfe')
        
        assert df.shape[0] == 10
        assert len(idx) == 10
        
    def test_correct_input(self):
        '''
        Test that only FAFSA data can pass through function.
        '''
        
        test_df = pd.DataFrame({'a':[1,2,3], 'b':[4,5,6]})
        
        with self.assertRaises(AttributeError):
            c2.fafsa_sample(test_df, size = 10)

        
class TestFrankWolfe(unittest.TestCase):
    ''' 
    Tests for the FrankWolfe function.
    '''
    
    def test_smoke(self):
        '''
        Smoke test to ensure function runs without error.
        '''
        assert c2.FrankWolfe(np.random.rand(3,100), 0.01)
        
    def test_oneshot(self):
        '''
        Test that reasonable probability distribution is returned from known inputs.
        '''
        lamb, _, _ = c2.FrankWolfe(np.random.rand(4,50), 0.01)
        assert np.isclose(np.sum(lamb), 1)

        
class TestgPrime(unittest.TestCase):
    ''' 
    Tests for the g_prime_i function.
    '''
    
    def test_smoke(self):
        '''
        Smoke test to ensure function runs without error.
        '''
        X = np.array([[1,2,3],[2,3,4],[3,4,5]])
        A = np.array([[1,2,3],[2,3,4],[3,4,5]])
        i = 0
        c2.g_prime_i(X, A, i)
        
    def test_oneshot(self):
        '''
        Test that function returns correct values with known input.
        '''
        X = np.array([[1,2,3],[2,3,4],[3,4,5]])
        A = np.array([[1,2,3],[2,3,4],[3,4,5]])
        i = 2
        assert c2.g_prime_i(X, A, i)[0][0] == 22.0      