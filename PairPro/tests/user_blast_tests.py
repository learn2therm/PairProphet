import unittest
import os
import pandas as pd

from PairPro.user_blast import make_blast_df

def get_test_df(filename):
    
    df_path = os.path.abspath(os.path.join('..', '..', 'data', filename))
    df = pd.read_csv(df_path)
    
    return df
    
class TestUserBlast(unittest.TestCase):
    '''
    Tests for the connect_db function.
    '''

    def test_smoke(self):
        '''
        Smoke test to ensure alignments can be made without error
        '''

        df = get_test_df('make_blast_df_test.csv')
        make_blast_df(df)


    def test_df(self):
        '''
        Test to make sure a DataFrame is the function output.
        '''
