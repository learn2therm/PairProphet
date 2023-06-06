"""
Unit tests for input cleaning and ML modules.
"""
import unittest
import duckdb as ddb
import pandas as pd
import numpy as np
import os

from pairpro.train_val_classification import train_model
from pairpro.train_val_classification import validate_model
from pairpro.train_val_classification import rf_wrapper

from pairpro.train_val_input_cleaning import check_input_type
from pairpro.train_val_input_cleaning import clean_input_columns
from pairpro.train_val_input_cleaning import verify_input_columns
from pairpro.train_val_input_cleaning import check_input_nans
from pairpro.train_val_input_cleaning import verify_protein_pairs
from pairpro.train_val_featuregen import get_fasta_from_dataframe
from pairpro.train_val_featuregen import get_protein_descriptors
from pairpro.train_val_featuregen import clean_new_dataframe
from pairpro.train_val_featuregen import create_new_dataframe


# conn = ddb.connect('pairpro/tests/l2t_mini_pdb.db', read_only=True)

# df = con.execute(f"""SELECT pair_id, m_protein_seq, t_protein_seq, bit_score, local_gap_compressed_percent_id, 
#         scaled_local_query_percent_id, scaled_local_symmetric_percent_id, 
#         query_align_len, query_align_cov, subject_align_len, subject_align_cov, 
#         LENGTH(m_protein_seq) AS m_protein_len, LENGTH(t_protein_seq) AS t_protein_len, FROM """).df()

# unit tests for input cleaning

class TestInputType(unittest.TestCase):
    """
    Tests that input data is a pandas dataframe
    """

    def setUp(self):
        self.conn = ddb.connect('pairpro/tests/l2t_mini_pdb.db', read_only=True)
        return

    def test_input_type(self):
        try:
            check_input_type([4, 3])
            self.assertTrue(False)
        except AssertionError:
            self.assertTrue(True)

    def tearDown(self):
        return super().tearDown()

class TestInputCleaning(unittest.TestCase):
    """
    Checks that data has been cleaned properly. Asserts unwanted
    columns have been removed and necessary columns remain in the
    dataframe.
    """
    # pass through some titles that should not be in the dataframe

    def test_input_cleaning(self):
        for title in ['Unnamed: 0', 'm_seq', 't_seq', 'prot_pair_index']:
            assert title not in clean_input_columns(df)

    # try to drop a column that should not be there
    def test_column_verification(self):
        try:
            verify_input_columns(df.drop(columns='meso_ogt'))
            self.assertTrue(False)
        except KeyError:
            self.assertTrue(True)


class TestForNans(unittest.TestCase):
    """
    Checks for NaN values in the dataframe and removes rows
    with NaN values.
    """

    def test_input_nans(self):
        df['another_column'] = pd.DataFrame([np.nan for i in range(len(df))])
        assert check_input_nans(df).isna().any().any() == False


class TestProteinPairs(unittest.TestCase):
    """
    Makes sure that there is a mesophillic and thermophyllic sequence for
    each row.
    """

    def test_protein_pair(self):
        try:
            verify_protein_pairs(
                df.drop(columns=(['m_protein_len', 't_protein_len'])))
            self.assertTrue(False)
        except AssertionError:
            self.assertTrue(True)

# unit tests for ML model


input_features = [
    'local_gap_compressed_percent_id',
    'scaled_local_query_percent_id',
    'scaled_local_symmetric_percent_id',
    'query_align_len',
    'query_align_cov',
    'subject_align_len',
    'subject_align_cov',
    'bit_score',
    'm_protein_len',
    't_protein_len']

target = 'protein_match'


class TestModelTraining(unittest.TestCase):
    """
    Checks input and output formats for model training
    function.
    """

    def test_invalid_inputs(self):
        # test that input data type is correct
        try:
            train_model([1, 2, 3], columns='string', target='string')
            self.assertTrue(False)
        except AssertionError:
            self.assertTrue(True)

    def test_output_format(self):
        """
        Asserts that function returns 5 objects to be assigned to
        model, train_X, train_y, dev_X, dev_y.
        """
        assert len(train_model(df, columns=input_features,
                               target=target)) == 5


class TestModelPerformance(unittest.TestCase):
    """
    Asserts that model testing function is taking in
    a pandas dataframe and outputting a numpy array.
    Also checks that the length of the predictions
    vector matches the dimensions of the target vector.
    """

    def test_asserts(self):
        model, _, _, val_X, val_y = train_model(
            df, columns=input_features, target='protein_match'
        )
        # assert that input types are correct
        with self.assertRaises(AssertionError):
            validate_model(model, [1, 2, 3], val_y)

    def test_model_output(self):
        model, _, _, val_X, val_y = train_model(
            df, columns=input_features, target='protein_match'
        )
        # assert output type is correct
        output, _ = validate_model(model, val_X, val_y)
        self.assertIsInstance(output, np.ndarray)

    def test_pred_dimension(self):
        model, _, _, val_X, val_y = train_model(
            df, columns=input_features, target='protein_match')
        # want to check that the number of predictions is equal to the number
        # of test examples
        preds, _ = validate_model(model, val_X, val_y)
        self.assertEqual(len(val_y), len(preds))


class TestWrapper(unittest.TestCase):
    """
    Ensures that the ML wrapping function is executing in the
    correct order.
    """

    def test_wrapper_input(self):
        # test that input data type is correct
        try:
            rf_wrapper([1, 2, 3], target=['hmmer_match'])
            self.assertTrue(False)
        except AssertionError:
            self.assertTrue(True)

    def test_wrapper_output(self):
        model, _, _, val_X, val_y = train_model(
            df,
            columns=input_features, target='protein_match'
        )
        # assert output type is correct
        output, _ = validate_model(model, val_X, val_y)
        self.assertIsInstance(output, np.ndarray)

    def test_output_dimension(self):
        model, _, _, val_X, val_y = train_model(
            df,
            columns=input_features, target='protein_match'
        )
        # want to check that the # of predictions is equal to # of examples
        preds, _ = validate_model(model, val_X, val_y)
        self.assertEqual(len(val_y), len(preds))


"""
Test feature generation functions
"""

class TestGetFastaFromDataframe(unittest.TestCase):

    def test_get_fasta_from_dataframe(self):
        # Example input dataframe
        dataframe = pd.DataFrame({
            'prot_pair_index': [1, 2, 3],
            'm_protein_seq': ['MESOSEQ1', 'MESOSEQ2', 'MESOSEQ3'],
            't_protein_seq': ['THERMOSEQ1', 'THERMOSEQ2', 'THERMOSEQ3']
        })

        # Set up output file paths
        output_file_a = 'output_a.fasta'
        output_file_b = 'output_b.fasta'

        try:
            # Call the function
            fasta_files = get_fasta_from_dataframe(dataframe, output_file_a, output_file_b)

            # Assert the existence of output files
            self.assertTrue(os.path.exists(output_file_a))
            self.assertTrue(os.path.exists(output_file_b))

            # Read the contents of the output files
            with open(output_file_a, 'r') as file_a:
                contents_a = file_a.read()
            with open(output_file_b, 'r') as file_b:
                contents_b = file_b.read()

            # Assert the contents of the output files
            expected_a = '>1\nMESOSEQ1\n>2\nMESOSEQ2\n>3\nMESOSEQ3\n'
            expected_b = '>1\nTHERMOSEQ1\n>2\nTHERMOSEQ2\n>3\nTHERMOSEQ3\n'
            self.assertEqual(contents_a, expected_a)
            self.assertEqual(contents_b, expected_b)

            # Assert the returned file paths
            self.assertEqual(fasta_files, [output_file_a, output_file_b])

        finally:
            # Clean up - delete the output files
            if os.path.exists(output_file_a):
                os.remove(output_file_a)
            if os.path.exists(output_file_b):
                os.remove(output_file_b)

