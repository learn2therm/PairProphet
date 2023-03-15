import pandas as pd
import unittest
from find_structures import scrapepdb

class ScrapePdbTest(unittest.TestCase):
    def test_return_polymer_entity(self):
        """
        Test if returned result from query is a Python dict. 
        If not, an AssertionError message is showed
        """
