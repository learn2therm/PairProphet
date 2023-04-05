"""
Unit tests for input cleaning and ML modules.
"""
import unittest
import pandas as pd
import numpy as np

from FAFSA.train_val_classification import train_model
from FAFSA.train_val_classification import evaluate_model
from FAFSA.train_val_classification import rf_wrapper

from FAFSA.train_val_input_cleaning import check_input_type
from FAFSA.train_val_input_cleaning import clean_input_columns
from FAFSA.train_val_input_cleaning import verify_input_columns
from FAFSA.train_val_input_cleaning import check_input_nans
from FAFSA.train_val_input_cleaning import verify_protein_pairs

df = pd.read_csv('learn2therm_sample_50k.csv')

# create new Boolean column for protein functionality match
df['protein_match'] = df['t_protein_desc'] == df['m_protein_desc']


# unit tests for input cleaning

class TestInputType(unittest.TestCase):
    """
    Tests that input data is a pandas dataframe
    """

    def test_input_type(self):
        try:
            check_input_type([4, 3])
            self.assertTrue(False)
        except AssertionError:
            self.assertTrue(True)


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
            evaluate_model(model, [1, 2, 3], val_y)

    def test_model_output(self):
        model, _, _, val_X, val_y = train_model(
            df, columns=input_features, target='protein_match'
        )
        # assert output type is correct
        output, _ = evaluate_model(model, val_X, val_y)
        self.assertIsInstance(output, np.ndarray)

    def test_pred_dimension(self):
        model, _, _, val_X, val_y = train_model(
            df, columns=input_features, target='protein_match')
        # want to check that the number of predictions is equal to the number
        # of test examples
        preds, _ = evaluate_model(model, val_X, val_y)
        self.assertEqual(len(val_y), len(preds))


class TestWrapper(unittest.TestCase):
    """
    Ensures that the ML wrapping function is executing in the
    correct order.
    """

    def test_wrapper_input(self):
        # test that input data type is correct
        try:
            rf_wrapper([1, 2, 3])
            self.assertTrue(False)
        except AssertionError:
            self.assertTrue(True)

    def test_wrapper_output(self):
        model, _, _, val_X, val_y = train_model(
            df,
            columns=input_features, target='protein_match'
        )
        # assert output type is correct
        output, _ = evaluate_model(model, val_X, val_y)
        self.assertIsInstance(output, np.ndarray)

    def test_output_dimension(self):
        model, _, _, val_X, val_y = train_model(
            df,
            columns=input_features, target='protein_match'
        )
        # want to check that the # of predictions is equal to # of examples
        preds, _ = evaluate_model(model, val_X, val_y)
        self.assertEqual(len(val_y), len(preds))
