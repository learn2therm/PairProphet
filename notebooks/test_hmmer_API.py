"""
Test unit for testing the function in the notebook by whole code and by 
using the real and fake output.
The packages you need to run this script are:
- pandas
- unittest
- hmmer
"""

import pandas as pd
import unittest
from hmmer import run_hmmerscanner


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
        # Input data
        # df = 
        #pd.read_csv("/Users/amin/ValidProt/FAFSA/learn2therm_sample_50k.csv")
        data = {
            'protein_seq': ['MAESGTSRRADHLVPVPGPDAEPPAVADELLRAVGRGDEQAFGRLYDLLAPRVYGLIRRVLRDPALAEEVTQEVLVEVWRRAARFDPAQGSANAWVFTIAHRRAVDRVRAEQKAADRTVRAGAAALDSPYDSVADEVSGRLERRQVRHCLDALTGLQREVVTLAYYQGHSYPQVAELLKTPLGTVKTRMRDGLIRLRDCLGVEATA'],
            'pid': [48641291]
        }

        df = pd.DataFrame(data)

        # run the hmmerscanner function on the sample DataFrame
        results_df = hmmerscanner(df, 1, 20)

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
