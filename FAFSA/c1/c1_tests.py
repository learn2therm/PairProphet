import pandas as pd
import numpy as np
import time

import unittest
import duckdb

import os
import sys

import c1

def get_db_path(filename = 'validprot_testing'):
    '''
    Gets path to unit test dataset for testing functions.

    Args:
        filename (str): Name of the unit test database file. Generally should not be
        altered. validprot_testing should allow complete and efficient unit tests for
        the full learn2therm database, but is only 1.5 MB and will run through tests
        very quickly on most systems.

    Returns:
        db_path (str): Full path to unit test database file.

    Raises:
        ValueError: filename not found in current directory.
    '''

    # Get path for test dataset import
    db_path = os.path.abspath(os.path.join('../', filename))

    if os.path.exists(db_path) is False:
        raise ValueError(f'Could not find {filename} in current directory')

    return db_path


class TestFetch(unittest.TestCase):
    ''' 
    Tests for the fetch_data function.
    '''
    
    def test_smoke(self):
        '''
        Smoke test to ensure DataFrame can be build without error.
        '''
        db_path = get_db_path()
        c1.fetch_data(db_path)
        
    def test_oneshot(self):
        '''
        Oneshot tests for supported combinations of file type and sampling method.
        '''
        assert False
        
        
    def test_wrong_type(self):
        '''
        Test for unsupported file type input.
        '''
        assert False
        
        
    def test_wrong_method(self):
        '''
        Test for unsupported sampling method input.
        '''
        assert False
        
    def test_idx_out_of_range(self):
        '''
        Test for unusable index range.
        '''
        assert False
        
    
        
        
    

