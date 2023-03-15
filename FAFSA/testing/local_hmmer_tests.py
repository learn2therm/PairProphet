""" 
tests for component three.

Note: the second tests that are commented out work, 
but because the nature of compute_local_hmmer.py changed, 
there are no outputs to check for
"""
# system dependecies
import glob
import os
import unittest

# library dependencies
import pandas as pd

# biopython
from Bio import SeqIO


# local dependencies/utils
from FAFSA.compute_local_hmmer import read_seq



class TestReadSeq(unittest.TestCase):
    """
    Tests for read_seq function
    """

    # inputs in setup function
    def setUp(self):
        self.sequences = pd.DataFrame(
            {
                'protein_seq': [
                    'MTITLVEGSEQLRKQVAYVEELRS',
                    'MFRDYVYDYLMRLWRVEYRH',
                    'MVTVLHESDQSLQGQAYLAEELRS',
                    'MTRTLVEGSEQLRKQVAYVEELRS']})
        self.inputname = 'test_input'
        self.empty_seq = pd.DataFrame({'protein_seq': []})
        self.invalid_seq = pd.DataFrame({'protein_seq': ['AAGCTYY']})
        read_seq(self.sequences, self.inputname)

        # create a new empty dataframe instance
        self.empty_df = pd.DataFrame()

    # Methods/tests

    def test_file_created(self):
        """ Tests if the file input is created"""
        self.assertTrue(os.path.isfile(f"{self.inputname}.fasta"))

    def test_file_contains_correct_sequences(self):
        """ Tests if the input is formatted correctly"""
        read_seq(self.sequences, self.inputname)
        records = list(SeqIO.parse(f"{self.inputname}.fasta", "fasta"))
        self.assertEqual(len(records), len(self.sequences))
        for i, record in enumerate(records):
            self.assertEqual(str(record.seq),
                             self.sequences.iloc[i]['protein_seq'])

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
        with self.assertRaises(ValueError) as err:
            read_seq(self.empty_df, self.inputname)
        self.assertEqual(str(err.exception), "Input dataframe is empty")

    def test_invalid_sequence(self):
        """Tests if the function raises an error when given an invalid sequence"""
        with self.assertRaises(AttributeError) as err:
            read_seq(self.invalid_seq, self.inputname)
            print(err.exception)

    # tear down of set up

    def tearDown(self):
        try:
            os.remove(f"{self.inputname}.fasta")
        except FileNotFoundError:
            pass


# class TestRunHmmer(unittest.TestCase):
#     """Tests for run_hmmer function"""
#     # set up

#     def setUp(self):
#         file_list = glob.glob("meso_output*.domtblout")
#         if len(file_list) > 0:
#             self.outputname = file_list[0][:-10]
#         else:
#             file_list = glob.glob("thermo_output*.domtblout")
#             if len(file_list) > 0:
#                 self.outputname = file_list[0][:-10]
#             else:
#                 self.outputname = None
#     # check if file is created. I don't need more tests because HMMER isn't my
#     # function

#     def test_file_created(self):
#         """Tests if the file output is created"""
#         if self.outputname:
#             self.assertTrue(os.path.isfile(f"{self.outputname}.domtblout"))
#         else:
#             self.fail(
#                 "No file matching meso_output*.domtblout or thermo_output*.domtblout was found")
