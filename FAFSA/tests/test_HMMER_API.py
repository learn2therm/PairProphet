"""
Test unit for testing the function in the notebook by whole code and by using the real and fake output.
The packages you need to run this script are:
- pandas
- unittest
- API
"""

import pandas as pd
import unittest
from HMMER_API import hmmerscanner


class HmmerScannerTest(unittest.TestCase):
    """
    This code defines a unit test class named `HmmerScannerTest` that inherits
    from `unittest.TestCase`. It includes two test methods:
    `test_realData_hmmerscanner` and `test_fakeData_hmmerscanner`.
    """

    def test_realData_hmmerscanner(self):
        """
        In `test_realData_hmmerscanner`, the hmmerscanner function is run on a
        sample DataFrame, and assertions are made to check if the output
        DataFrame contains the expected columns. If an assertion fails,
        an error message is printed to indicate which assertion failed.
        """
        # Read the input data and print here
        df = pd.read_csv("/Users/amin/ValidProt/FAFSA/learn2therm_sample_50k.csv")

        # run the hmmerscanner function on the sample DataFrame
        results_df = hmmerscanner(df, 1)

        # assertion 1: check if the output DataFrame is not empty
        try:
            assert len(results_df) > 0
            print("assertion 1 passed")
        except AssertionError:
            print("assertion 1 failed, because the HMMER server may have a problem")

        # assertion 2: check if the output DataFrame contains 'acc' column
        try:
            assert 'acc' in results_df.columns
            print("assertion 2 passed")
        except AssertionError:
            print("assertion 2 failed, because the HMMER server may have a problem")

        # assertion 3: check if the output DataFrame contains 'name' column
        try:
            assert 'name' in results_df.columns
            print("assertion 3 passed")
        except AssertionError:
            print("assertion 3 failed, because the HMMER server may have a problem")

        # assertion 4: check if the output DataFrame contains 'score' column
        try:
            assert 'score' in results_df.columns
            print("assertion 4 passed")
        except AssertionError:
            print("assertion 4 failed, because the HMMER server may have a problem")

        # assertion 5: check if the output DataFrame contains 'evalue' column
        try:
            assert 'evalue' in results_df.columns
            print("assertion 5 passed")
        except AssertionError:
            print("assertion 5 failed, because the HMMER server may have a problem")

        # assertion 6: check if the output DataFrame contains 'pvalue' column
        try:
            assert 'pvalue' in results_df.columns
            print("assertion 6 passed")
        except AssertionError:
            print("assertion 6 failed, because the HMMER server may have a problem")

        # assertion 7: check if the output DataFrame contains 'desc' column
        try:
            assert 'desc' in results_df.columns
            print("assertion 7 passed")
        except AssertionError:
            print("assertion 7 failed, because the HMMER server may have a problem")

        # assertion 8: check if the output DataFrame does not contain 'tlen' column
        try:
            assert 'tlen' not in results_df.columns
            print("assertion 8 passed")
        except AssertionError:
            print("assertion 8 failed, because the HMMER server may have a problem")

        # assertion 9: check if the output DataFrame does not contain 'ali_len' column
        try:
            assert 'ali_len' not in results_df.columns
            print("assertion 9 passed")
        except AssertionError:
            print("assertion 9 failed, because the HMMER server may have a problem")

        # assertion 10: check if the output DataFrame does not contain 'env_from' column
        try:
            assert 'env_from' not in results_df.columns
            print("assertion 10 passed")
        except AssertionError:
            print("assertion 10 failed, because the HMMER server may have a problem")

    def test_fakeData_hmmerscanner(self):
        """
        In `test_fakeData_hmmerscanner`, a fake DataFrame is created, and assertions are made to
        check if it has the expected columns. If an assertion fails, an error message is printed.
        It is with fake data for testing the code when we get the error because of the server problem.
        """
        # create a fake DataFrame
        results_dff = pd.DataFrame({
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

        # assertion 1: check if the output is a DataFrame
        try:
            assert isinstance(results_dff, pd.DataFrame)
            print("assertion 1 passed")
        except AssertionError:
            print("assertion 1 failed")

        # assertion 2: check if the output DataFrame has the expected columns
        try:
            expected_columns = ['sequence', 'hmm_acc', 'hmm_desc', 'hmm_from', 'hmm_to', 'hmm_cov', 'hmm_score',
                                'hmm_evalue', 'hmm_pvalue', 'acc', 'name', 'score', 'evalue', 'pvalue', 'desc']
            assert set(results_dff.columns) == set(
                expected_columns), "Unexpected columns in output DataFrame"
            print("assertion 2 passed")
        except AssertionError:
            print("assertion 2 failed")
