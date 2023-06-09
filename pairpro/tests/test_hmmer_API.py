"""
Test unit for testing the function in the notebook by whole code and by 
using the real and fake output.
The packages you need to run this script are:
- pandas
- unittest
- hmmer
"""

import pandas as pd
import asyncio
import unittest
from pairpro.hmmer import run_hmmerscanner

class HmmerScannerTest(unittest.TestCase):
    """
    This code defines a unit test class named `HmmerScannerTest` that inherits
    from `unittest.TestCase`. It includes a test method:
    `test_realData_hmmerscanner`.
    """

    def test_realData_hmmerscanner(self):
        """
        In `test_realData_hmmerscanner`, the hmmerscanner function is run on a
        sample DataFrame, and assertions are made to check if the output
        DataFrame contains the expected columns. If an assertion fails,
        an error message is printed to indicate which assertion failed.
        """

        # Input data
        data = {
            'protein_seq': ['MAESGTSRRADHLVPVPGPDAEPPAVADELLRAVGRGDEQAFGRLYDLLAPRVYGLIRRVLRDPALAEEVTQEVLVEVWRRAARFDPAQGSANAWVFTIAHRRAVDRVRAEQKAADRTVRAGAAALDSPYDSVADEVSGRLERRQVRHCLDALTGLQREVVTLAYYQGHSYPQVAELLKTPLGTVKTRMRDGLIRLRDCLGVEATA'],
            'pair_id': [48641291]
        }

        df = pd.DataFrame(data)

        # run the hmmerscanner function on the sample DataFrame
        results_df = run_hmmerscanner(df, "protein_seq", 1, 20, '')

        # assertion 1: check if the output DataFrame is not empty
        assert len(
            results_df) > 0, "assertion 1 failed, because the HMMER server may have a problem"

        # assertion 2: check if the output DataFrame contains 'acc' column
        assert 'acc' in results_df.columns, "assertion 2 failed, because the HMMER server may have a problem"

        # assertion 3: check if the output DataFrame contains 'name' column
        assert 'name' in results_df.columns, "assertion 3 failed, because the HMMER server may have a problem"

        # assertion 4: check if the output DataFrame contains 'score' column
        assert 'score' in results_df.columns, "assertion 4 failed, because the HMMER server may have a problem"

        # assertion 5: check if the output DataFrame contains 'evalue' column
        assert 'evalue' in results_df.columns, "assertion 5 failed, because the HMMER server may have a problem"

        # assertion 6: check if the output DataFrame contains 'pvalue' column
        assert 'pvalue' in results_df.columns, "assertion 6 failed, because the HMMER server may have a problem"

        # assertion 7: check if the output DataFrame contains 'desc' column
        assert 'desc' in results_df.columns, "assertion 7 failed, because the HMMER server may have a problem"

        # assertion 8: check if the output DataFrame does not contain 'tlen' column
        assert 'tlen' not in results_df.columns, "assertion 8 failed, because the HMMER server may have a problem"

        # assertion 9: check if the output DataFrame does not contain 'ali_len' column
        assert 'ali_len' not in results_df.columns, "assertion 9 failed, because the HMMER server may have a problem"

        # assertion 10: check if the output DataFrame does not contain 'env_from' column
        assert 'env_from' not in results_df.columns, "assertion 10 failed, because the HMMER server may have a problem"


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(HmmerScannerTest)
    unittest.TextTestRunner().run(suite)
