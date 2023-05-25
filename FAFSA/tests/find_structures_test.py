"""
This module contains test for find_structures.py
"""

import pandas as pd
import unittest
from FAFSA.find_structures import scrapepdb

seq = "MLLSDRDLVSEIKSGDLSLEPFEPALLQPSSIDVRLDRFFRVFNNHLYTHIDPAEQQDDLTAEVEVTDGEAFVLHPGEFVLASTLEVITLGDQLAGRLEGKSSLGRLGLLTHSTAGFIDPGFSGHVTLELSNVANLPIKLWPGMKIGQLCIFRLSSPAEHPYGSAVYGSRYQGQRGPTPSRSAQNFRLWPTS"
df_result_reind = scrapepdb(seq)

class ScrapePdbTest(unittest.TestCase):
    """
    Tests for scrapepdb()
    """
    def test_invalid_output(self):
        assert type(df_result_reind) == pd.DataFrame, f"Unexpected type returned"
        assert len(df_result_reind) > 0, f"Invalid return"
        assert len(df_result_reind.columns) == 5, f"Expected 5 columns after formating, but got {len(df_result_reind.columns)} columns"