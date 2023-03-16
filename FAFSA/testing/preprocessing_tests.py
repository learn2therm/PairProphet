'''
Unit testing script for FAFSA preprocessing.
'''

import unittest

import os

from FAFSA.preprocessing import connect_db
from FAFSA.preprocessing import build_fafsa
from FAFSA.preprocessing import sankey_plots

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
    db_path = os.path.abspath(os.path.join('..', '..', 'data', filename))

    if os.path.exists(db_path) is False:
        raise ValueError(f'Could not find {filename} in current directory')

    return db_path



class TestConnection(unittest.TestCase):
    '''
    Tests for the connect_db function.
    '''

    def test_smoke(self):
        '''
        Smoke test to ensure connection can be made to database without error
        '''

        db_path = get_db_path()
        connect_db(db_path)


    def test_empty(self):
        '''
        Test to make sure an empty database triggers an AssertionError.
        '''

        db_path = get_db_path('fafsa_empty')

        with self.assertRaises(AttributeError):
            connect_db(db_path)


class TestBuildFafsa(unittest.TestCase):
    '''
    Tests for the build_fafsa function.
    '''

    def test_smoke(self):
        '''
        Smoke test to ensure database is assembled without error.
        '''

        db_path = get_db_path()
        con = connect_db(db_path)

        build_fafsa(con)

    def test_oneshot(self):
        '''
        One shot test with known FAFSA test set.
        '''

        db_path = get_db_path()
        con = connect_db(db_path)

        build_fafsa(con)
        con = connect_db(db_path)

        correct = ['proteins', 'protein_pairs', 'taxa', 'taxa_pairs', 'vp_final',
                   'vp_proteins', 'vp_protein_pairs', 'vp_taxa', 'vp_taxa_pairs', 
                   'vp_ogt_taxa_pairs']

        sizes = []
        for c in correct:
            sizes.append(con.execute(f"""SELECT COUNT(*) FROM {c}""").df().values.flatten()[0])

        correct_sizes = [400, 200, 250, 199, 269, 400, 263, 250, 199, 253]

        assert sizes == correct_sizes

    def test_db_format(self):
        '''
        Tests that genuine learn2therm input generates correct FAFSA tables.
        '''

        db_path = get_db_path()
        con = connect_db(db_path)

        tables = con.execute("""SELECT TABLE_NAME
                                FROM INFORMATION_SCHEMA.TABLES
                                WHERE TABLE_TYPE='BASE TABLE'""").df()
        correct = ['proteins', 'protein_pairs', 'taxa', 'taxa_pairs', 'vp_final',
                   'vp_proteins', 'vp_protein_pairs', 'vp_taxa', 'vp_taxa_pairs', 
                   'vp_ogt_taxa_pairs']
        assert (item in tables for item in correct)

    def test_ogt_out_of_range(self):
        '''
        Test for improper ogt spread.
        '''

        db_path = get_db_path()
        con = connect_db(db_path)

        with self.assertRaises(ValueError):
            build_fafsa(con, min_ogt_diff = -10)

    def test_16s_out_of_range(self):
        '''
        Test for ValueError on 16S cutoff below 1 bp.
        '''

        db_path = get_db_path()
        con = connect_db(db_path)

        with self.assertRaises(ValueError):
            build_fafsa(con, min_16s = 0)


class TestSankey(unittest.TestCase):
    '''
    Tests for sankey_plots function.
    '''

    def test_smoke(self):
        '''
        Smoke test to ensure Sankey plots are made without error.
        '''

        db_path = get_db_path()
        con = connect_db(db_path)

        min_ogt_diff = 20

        sankey_plots(con, min_ogt_diff)

    def test_ogt_out_of_range(self):
        '''
        Test for improper ogt spread.
        '''
        db_path = get_db_path()
        con = connect_db(db_path)

        min_ogt_diff = -10

        with self.assertRaises(ValueError):
            sankey_plots(con, min_ogt_diff)
