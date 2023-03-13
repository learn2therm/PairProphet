import unittest
import pandas as pd
from C3API import hmmerscanner


class HmmerScannerTest(unittest.TestCase):
    
    def test_realData_hmmerscanner1(self):
        
        
        # read CSV file
        df = pd.read_csv("/Users/amin/ValidProt/data/Sample.csv")
        

        # run the hmmerscanner function on the sample DataFrame
        results_df = hmmerscanner(df)

        # run assertions on the output DataFrame
        try:
            assert len(results_df) > 0
            print("assertion 1 passed")
        except AssertionError:
            print("assertion 1 failed")

        try:
            assert 'acc' in results_df.columns
            print("assertion 2 passed")
        except AssertionError:
            print("assertion 2 failed")

        try:
            assert 'name' in results_df.columns
            print("assertion 3 passed")
        except AssertionError:
            print("assertion 3 failed")

        try:
            assert 'score' in results_df.columns
            print("assertion 4 passed")
        except AssertionError:
            print("assertion 4 failed")

        try:
            assert 'evalue' in results_df.columns
            print("assertion 5 passed")
        except AssertionError:
            print("assertion 5 failed")

        try:
            assert 'pvalue' in results_df.columns
            print("assertion 6 passed")
        except AssertionError:
            print("assertion 6 failed")

        try:
            assert 'desc' in results_df.columns
            print("assertion 7 passed")
        except AssertionError:
            print("assertion 7 failed")

        try:
            assert 'tlen' not in results_df.columns
            print("assertion 8 passed")
        except AssertionError:
            print("assertion 8 failed")

        try:
            assert 'ali_len' not in results_df.columns
            print("assertion 9 passed")
        except AssertionError:
            print("assertion 9 failed")

        try:
            assert 'env_from' not in results_df.columns
            print("assertion 10 passed")
        except AssertionError:
            print("assertion 10 failed")
            
    def test_fakeData_hmmerscanner2(self):
        # It is with fake data for testing the code when we get the error because of the server problem
        # create a fake DataFrame
        results_df = pd.DataFrame({
            'sequence': ['AAAAAA', 'AAAAAA', 'BBBBBB'],
            'hmm_acc': ['PF00001', 'PF00002', 'PF00001'],
            'hmm_desc': ['Test domain 1', 'Test domain 2', 'Test domain 1'],
            'hmm_from': [1, 10, 20],
            'hmm_to': [50, 60, 70],
            'hmm_cov': [0.5, 0.6, 0.7],
            'hmm_score': [50.0, 60.0, 70.0],
            'hmm_evalue': [1e-6, 1e-7, 1e-8],
            'hmm_pvalue': [1e-8, 1e-9, 1e-10],
            'acc': ['ABC123', 'DEF456', 'GHI789'],
            'name': ['Protein 1', 'Protein 2', 'Protein 3'],
            'score': [100, 200, 300],
            'evalue': [1e-10, 1e-11, 1e-12],
            'pvalue': [1e-12, 1e-13, 1e-14],
            'desc': ['Protein description 1', 'Protein description 2', 'Protein description 3']
        })

        # test the output
        try:
            assert isinstance(results_df, pd.DataFrame)
            print("assertion 1 passed")
        except AssertionError:
            print("assertion 1 failed")
    
        try:
            assert set(results_df.columns) == set(['sequence', 'hmm_acc', 'hmm_desc', 'hmm_from', 'hmm_to', 'hmm_cov',                                           'hmm_score', 'hmm_evalue', 'hmm_pvalue', 'acc', 'name', 'score', 'evalue', 'pvalue', 'desc'])
            print("assertion 2 passed")
        except AssertionError:
            print("assertion 2 failed")
            

