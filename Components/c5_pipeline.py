import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.stats
import sklearn.preprocessing
import sklearn.model_selection
import sklearn.neighbors
import sklearn.ensemble
import sklearn.feature_selection

import c5_classification
import c5_input_cleaning

from c5_classification import train_model
from c5_classification import test_model
from c5_classification import plot_model
from c5_classification import rf_wrapper

from c5_input_cleaning import check_input_type
from c5_input_cleaning import clean_input_columns
from c5_input_cleaning import verify_input_columns
from c5_input_cleaning import check_input_NANs
from c5_input_cleaning import verify_protein_pairs
from c5_input_cleaning import input_cleaning_wrapper

# autopep8 --in-place --aggressive --aggressive c5_pipeline.py

df_original = pd.read_csv('learn2therm_sample_50k.csv')

"""
This cleaning step is done because we want to use thermophyllic protein
description as our target. I am eliminating categories without significant
value counts in order to build a better model with our "dummy" target.
This code will be deleted when we get our real target from Component 3.
"""

categories = df_original['t_protein_desc'].value_counts()
categories = categories.iloc[categories.values > 500]
categories_dict = {item: None for item in categories.index}
list_of_cats = list(categories_dict.keys())
list_of_cats

# process dataframe
df = df_original[df_original.t_protein_desc.isin(list_of_cats)]

df['protein_match'] = df['t_protein_desc'].eq(df['m_protein_desc'])

# will need to be adjusted when component 3 is finished
columns_to_keep = [
    'bit_score',
    'local_gap_compressed_percent_id',
    'scaled_local_query_percent_id',
    'scaled_local_symmetric_percent_id',
    'query_align_len',
    'query_align_cov',
    'subject_align_len',
    'subject_align_cov',
    'm_protein_len',
    't_protein_len',
    'protein_match']


# run dataframe through component five wrapper
def c5_wrapper(dataframe):

    # clean input dataframe
    cleaned = input_cleaning_wrapper(dataframe)

    # run through model
    model = rf_wrapper(cleaned)

    return model
