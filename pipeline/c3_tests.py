""" tests for component three"""
# system dependecies
from io import StringIO
import os
from pathlib import Path
import tempfile
from tempfile import TemporaryDirectory
import subprocess
import unittest

# library dependencies
from collections import defaultdict
import pandas as pd

## biopython
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SearchIO

# local dependencies/utils
from c3_pfam_hmmer import read_seq

## Paths
PFAM_PATH = Path("/Users/humoodalanzi/pfam/Pfam-A.hmm")
ID_DB_PATH = Path("/Users/humoodalanzi/pfam/proteins_id.zip")



class TestReadSeq(unittest.TestCase):
    """
    Tests for read_seq function
    """

    # inputs in setup function
    def setUp(self):
        self.sequences = pd.DataFrame({'protein_seq': ['MTITLVEGSEQLRKQVAYVEELRS',
                                           'MFRDYVYDYLMRLWRVEYRH',
                                           'MVTVLHESDQSLQGQAYLAEELRS',
                                           'MTRTLVEGSEQLRKQVAYVEELRS']})
        self.inputname = 'test_input'
        self.empty_seq = pd.DataFrame({'protein_seq': []})
        self.invalid_seq = pd.DataFrame({'protein_seq': ['AAGCTYY']})
        read_seq(self.sequences, self.inputname)

        # create a new empty dataframe instance
        self.empty_df = pd.DataFrame()
        
    ## Methods/tests
        
    def test_file_created(self):
        """ Tests if the file input is created"""
        self.assertTrue(os.path.isfile(f"{self.inputname}.fasta"))
        
    def test_file_contains_correct_sequences(self):
        """ Tests if the input is formatted correctly"""
        read_seq(self.sequences, self.inputname)
        records = list(SeqIO.parse(f"{self.inputname}.fasta", "fasta"))
        self.assertEqual(len(records), len(self.sequences))
        for i, record in enumerate(records):
            self.assertEqual(str(record.seq), self.sequences.iloc[i]['protein_seq'])
    
    def test_invalid_input(self):
        """Tests if the function raises an error with invalid input"""
        with self.assertRaises(ValueError):
            read_seq(self.empty_df, self.inputname)
    
    def test_invalid_type(self):
        """Tests if the function raises an error when given an invalid type"""
        with self.assertRaises(AttributeError):
            read_seq(123, self.inputname)
        
    def test_empty_sequences(self):
        """Tests if the function raises an error when given an empty DataFrame"""
        with self.assertRaises(ValueError) as cm:
            read_seq(self.empty_df, self.inputname)
        self.assertEqual(str(cm.exception), "Input dataframe is empty")
        
    def test_invalid_sequence(self):
        """Tests if the function raises an error when given an invalid sequence"""
        with self.assertRaises(AttributeError) as cm:
            read_seq(self.invalid_seq, self.inputname)
            print(cm.exception)


    # tear down of set up
    def tearDown(self):
        try:
            os.remove(f"{self.inputname}.fasta")
        except FileNotFoundError:
            pass
    