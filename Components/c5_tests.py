import unittest
import pandas as pd
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import seaborn as sns
import sklearn.preprocessing
import sklearn.model_selection
import sklearn.neighbors

import c5_classification
import c5_input_cleaning

from c5_classification import train_model
from c5_classification import test_model
from c5_classification import plot_model
from c5_classification import RF_wrapper

from c5_input_cleaning import check_input_type
from c5_input_cleaning import clean_input_columns
from c5_input_cleaning import verify_input_columns
from c5_input_cleaning import check_input_NANs
from c5_input_cleaning import verify_protein_pairs
from c5_input_cleaning import input_cleaning_wrapper

#import training dataframe to run all of the tests
df = pd.read_csv('learn2therm_sample_50k.csv')

"""
Unit tests for C5 input and data cleaning functions
"""

class TestInputType(unittest.TestCase):
    def test_input_type(self): 
        """
        Tests that input data is a pandas dataframe.
        """
        try:
            check_input_type([4,3])
            self.assertTrue(False)
        except AssertionError:
            self.assertTrue(True)

class TestInputCleaning(unittest.TestCase):
   #pass through some titles that should not be in the dataframe
    def test_input_cleaning(self):
        for title in ['Unnamed: 0','m_seq', 't_seq', 'prot_pair_index']:
            assert title not in clean_input_columns(df)
    #try to drop a column that should not be there
    
    def test_column_verification(self):
        try:
            verify_input_columns(df.drop(columns='meso_ogt'))
            self.assertTrue(False)
        except KeyError:
            self.assertTrue(True)

class TestForNans(unittest.TestCase):
    
    def test_input_Nans(self):
        df['another_column'] = pd.DataFrame([np.nan for i in range(len(df))])
        assert check_input_NANs(df).isna().any().any() == False

class TestProteinPairs(unittest.TestCase):
   
    def test_protein_pair(self):
        try:
            verify_protein_pairs(df.drop(columns = (['m_protein_len', 't_protein_len'])))
            self.assertTrue(False)
        except AssertionError:
            self.assertTrue(True)

"""
Unit tests for kNN model functions.
"""

input_features = ['local_gap_compressed_percent_id','scaled_local_query_percent_id',
 'scaled_local_symmetric_percent_id', 'query_align_len', 'query_align_cov',
 'subject_align_len', 'subject_align_cov', 'bit_score',
 'm_protein_len', 't_protein_len']

target = 'protein_match'

class TestModelTraining(unittest.TestCase):
    
    def test_invalid_inputs(self):
    
        #test that input data type is correct
        
        try:
            train_model([1,2,3], columns = 'string', target = 'string')
            self.assertTrue(False)
        except AssertionError:
            self.assertTrue(True)
        
    def test_output_format(self):
        
        #asserts that function returns 4 objects to be assigned to pearson_corr, model, test_X, test_y
        assert len(train_model(df, columns = input_features, 
                                          target = target)) == 5

class TestModelPerformance(unittest.TestCase):

    def test_asserts(self):
        model, dev_X, dev_y, test_X, test_y = train_model(
            df, columns=input_features, target='protein_match'
        )
        # assert that input types are correct
        with self.assertRaises(AssertionError):
            test_model(model, [1, 2, 3], test_y)
            
    def test_model_output(self):
        model, dev_X, dev_y, test_X, test_y = train_model(
            df, columns=input_features, target='protein_match'
        )
        # assert output type is correct
        output = test_model(model, test_X, test_y)
        self.assertIsInstance(output, np.ndarray)
        
    def test_pred_dimension(self):
        model, dev_X, dev_y, test_X, test_y = train_model(df, 
            columns=input_features, target='protein_match'
        )
        # want to check that the number of predictions is equal to the number of test examples
        preds = test_model(model, test_X, test_y)
        self.assertEqual(len(test_y), len(preds))

class TestWrapper(unittest.TestCase):
    
    def test_wrapper_input(self):
        #test that input data type is correct
        try:
            RF_wrapper([1,2,3])
            self.assertTrue(False)
        except AssertionError:
            self.assertTrue(True)

    def test_wrapper_output(self):
        model, dev_X, dev_y, test_X, test_y = train_model(
            df, 
            columns=input_features, target='protein_match'
        )
        # assert output type is correct
        output = test_model(model, test_X, test_y)
        self.assertIsInstance(output, np.ndarray)
        
    def test_output_dimension(self):
        model, dev_X, dev_y, test_X, test_y = train_model(
            df, 
            columns=input_features, target='protein_match'
        )
        # want to check that the number of predictions is equal to the number of test examples
        preds = test_model(model, test_X, test_y)
        self.assertEqual(len(test_y), len(preds))