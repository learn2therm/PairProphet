"""
To do: 
Write test for global vs local alignments
"""
import unittest
import os
import pandas as pd
import numpy as np
from PairPro.user_blast import make_blast_df

def get_test_df(filename):
    
    df_path = os.path.abspath(os.path.join('..', '..', 'data', filename))
    df = pd.read_csv(df_path)
    
    return df
    
class TestUserBlast(unittest.TestCase):
    '''
    Tests for the make_blast_df function.
    '''
    def test_smoke(self):
        '''
        Smoke test to ensure alignments can be made without error
        '''
        df = get_test_df('make_blast_df_test.csv')
        make_blast_df(df)

    def test_datatype(self):
        '''
        Test to make sure a DataFrame is the function output.
        '''
        df = get_test_df('make_blast_df_test.csv')
        blast_df = make_blast_df(df)
        self.assertTrue(type(blast_df) is pd.core.frame.DataFrame)
        
    def test_output_shape(self):
        '''
        Test to make sure output data has correct shape.
        '''
        df = get_test_df('make_blast_df_test.csv')
        blast_df = make_blast_df(df)
        self.assertTrue(blast_df.shape == (5,15))
        
    def test_col_names(self):
        '''
        Test to make sure output data has correct column names.
        '''
        cols = ['query', 'subject', 'pid', 'local_gap_compressed_percent_id',
                'scaled_local_query_percent_id', 'scaled_local_symmetric_percent_id',
                'query_align_len', 'query_align_cov', 'subject_align_len',
                'subject_align_cov', 'bit_score', 'query_id', 'subject_id']
        
        df = get_test_df('make_blast_df_test.csv')
        blast_df = make_blast_df(df)
        
        for col in cols:
            self.assertTrue(col in blast_df.columns)
            
    def test_metrics(self):
        '''
        Test to make sure metrics are calculated correctly.
        '''        
        df = get_test_df('make_blast_df_test.csv')
        blast_df = make_blast_df(df)
        
        correct_lgcpid = np.array([0.33472803, 0.72307692, 0.25555556, 0.27570093, 0.37083333])
        correct_slqpid = np.array([0.34934498, 0.72307692, 0.2754491 , 0.27830189, 0.3647541]) 
        correct_slspid = np.array([0.33898305, 0.72307692, 0.26512968, 0.27830189, 0.37316562])
        
        self.assertTrue(np.allclose(correct_lgcpid, blast_df['local_gap_compressed_percent_id'], atol=1e-04))
        self.assertTrue(np.allclose(correct_slqpid, blast_df['scaled_local_query_percent_id'], atol=1e-04))
        self.assertTrue(np.allclose(correct_slspid, blast_df['scaled_local_symmetric_percent_id'], atol=1e-04))
        
    def test_filter(self):
        '''
        Test to make sure improper rows are removed.
        '''        
        df = get_test_df('make_blast_df_test_empty_cell.csv')
        blast_df = make_blast_df(df)
       
        self.assertTrue(blast_df.shape == (3,15))