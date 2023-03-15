'''
Unit testing script for FAFSA c1.
'''

import unittest
import os

from FAFSA.model_input import *

def get_db_path(filename = 'fafsa_testing'):
    '''
    Gets path to unit test dataset for testing functions.

    Args:
        filename (str): Name of the unit test database file. Generally should not be
        altered. fafsa_testing should allow complete and efficient unit tests for
        the full learn2therm database, but is only 1.5 MB and will run through tests
        very quickly on most systems.

    Returns:
        db_path (str): Full path to unit test database file.

    Raises:
        ValueError: filename not found in current directory.
    '''

    # Get path for test dataset import
    db_path = os.path.abspath(os.path.join('..', 'data', filename))

    if os.path.exists(db_path) is False:
        raise ValueError(f'Could not find {filename} in current directory')

    return db_path


class TestFetch(unittest.TestCase):
    '''
    Tests for the fetch_data function.
    '''

    def test_smoke(self):
        '''
        Smoke test to ensure function runs in both csv and duckdb modes without error.
        '''
        db_path = get_db_path()
        c1.fetch_data(db_path, form = 'duckdb')
        df_path = get_db_path('learn2therm_sample_50k.zip')
        c1.fetch_data(df_path, form = 'csv', size = 100)

    def test_oneshot_random(self):
        '''
        Oneshot test for random sampling.
        '''
        df_path = get_db_path('learn2therm_sample_50k.zip')
        assert c1.fetch_data(df_path, form = 'csv', size = 100, method = 'random').shape[0] == 100

    def test_oneshot_chunk(self):
        '''
        Oneshot test for chunked sampling.
        '''
        df_path = get_db_path('learn2therm_sample_50k.zip')
        assert c1.fetch_data(df_path, form = 'csv', size = 100, method = 'chunk',
                      chunksize = 9).shape[0] == 100

    def test_oneshot_num(self):
        '''
        Oneshot test for numeric sampling. Test that size argument is superceded by idx_range.
        '''
        df_path = get_db_path('learn2therm_sample_50k.zip')
        assert c1.fetch_data(df_path, form = 'csv', size = 1, method = 'numeric',
                      idx_range = [100, 200]).shape[0] == 100

    def test_wrong_type(self):
        '''
        Test for unsupported file type input.
        '''
        db_path = get_db_path()

        with self.assertRaises(ValueError):
            c1.fetch_data(db_path, form = 'wrong')

    def test_wrong_method(self):
        '''
        Test for unsupported sampling method input.
        '''
        db_path = get_db_path()

        with self.assertRaises(ValueError):
            c1.fetch_data(db_path, method = 'wrong')

    def test_idx_out_of_range(self):
        '''
        Test for unusable index range.
        '''
        db_path = get_db_path()

        with self.assertRaises(ValueError):
            c1.fetch_data(db_path, method = 'sequential', idx_range = [0,10000])
