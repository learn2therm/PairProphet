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

from c5_classification import split_data
from c5_classification import train_reg
from c5_classification import test_reg
from c5_classification import plot_regression
from c5_classification import kNN_wrapper
from c5_classification import JSD_dev_and_test

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

class TestDataSplit(unittest.TestCase):
    
    def test_data_split_input(self):
        #test that input data type is correct
        try:
            split_data([1,2,3])
            self.assertTrue(False)
        except ValueError:
            self.assertTrue(True)
    
    def test_data_split_output(self):
        #test that output data type is correct
        if "tuple" in str(type(split_data(df))):
            self.assertTrue(True)
        else:
            self.assertTrue(False)

input_features = ['local_gap_compressed_percent_id','scaled_local_query_percent_id',
 'scaled_local_symmetric_percent_id', 'query_align_len', 'query_align_cov',
 'subject_align_len', 'subject_align_cov', 'bit_score',
 'm_protein_len', 't_protein_len']

target = 't_protein_desc'

class TestModelTraining(unittest.TestCase):
    
    def test_invalid_inputs(self):
        #test that input data type is correct
        try:
            train_reg([1,2,3], [4,5,6], columns = 'string', target = 'string')
            self.assertTrue(False)
        except AssertionError:
            self.assertTrue(True)
    
    def test_input_distro(self):
        #test that dev and test features have similar Jensen Shannon Divergence
        JSD = (
            scipy.stats.bootstrap((train_reg(split_data(df)[0], split_data(df)[1],
                                columns = input_features, target=target)[1], train_reg(split_data(df)[0], 
                                split_data(df)[1], columns = input_features, target=target)[3]), 
                                  JSD_dev_and_test, n_resamples=1000, 
                                  batch=5, method='percentile')
        )

        div = JSD.confidence_interval[1]
        
        #asserts that the divergence between data sets is sufficiently low
        assert abs(div) < 0.3, "Warning! High JSD between dev and test set!"
        
    def test_output_format(self):
        #asserts that function returns 4 objects to be assigned to pearson_corr, model, test_X, test_y
        assert len(train_reg(split_data(df)[0], split_data(df)[1], columns = input_features, 
                                          target = target)) == 5

class TestModelPerformance(unittest.TestCase):

    def test_asserts(self):
        model, dev_X, dev_y, test_X, test_y = train_reg(
            split_data(df)[0], split_data(df)[1], 
            columns=input_features, target='t_protein_desc'
        )
        # assert that input types are correct
        with self.assertRaises(AssertionError):
            test_reg(model, [1, 2, 3], test_y)
            
    def test_model_output(self):
        model, dev_X, dev_y, test_X, test_y = train_reg(
            split_data(df)[0], split_data(df)[1], 
            columns=input_features, target='t_protein_desc'
        )
        # assert output type is correct
        output = test_reg(model, test_X, test_y)
        self.assertIsInstance(output, np.ndarray)
        
    def test_pred_dimension(self):
        model, dev_X, dev_y, test_X, test_y = train_reg(
            split_data(df)[0], split_data(df)[1], 
            columns=input_features, target='t_protein_desc'
        )
        # want to check that the number of predictions is equal to the number of test examples
        preds = test_reg(model, test_X, test_y)
        self.assertEqual(len(test_y), len(preds))

class TestWrapper(unittest.TestCase):
    
    def test_wrapper_input(self):
        #test that input data type is correct
        try:
            kNN_wrapper([1,2,3])
            self.assertTrue(False)
        except AssertionError:
            self.assertTrue(True)

    def test_wrapper_output(self):
        model, dev_X, dev_y, test_X, test_y = train_reg(
            split_data(df)[0], split_data(df)[1], 
            columns=input_features, target='t_protein_desc'
        )
        # assert output type is correct
        output = test_reg(model, test_X, test_y)
        self.assertIsInstance(output, np.ndarray)
        
    def test_output_dimension(self):
        model, dev_X, dev_y, test_X, test_y = train_reg(
            split_data(df)[0], split_data(df)[1], 
            columns=input_features, target='t_protein_desc'
        )
        # want to check that the number of predictions is equal to the number of test examples
        preds = test_reg(model, test_X, test_y)
        self.assertEqual(len(test_y), len(preds))

