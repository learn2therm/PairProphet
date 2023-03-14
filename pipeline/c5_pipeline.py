"""
Wrapper functions for all of the machine learning component.
"""

import pandas as pd
from c5_classification import rf_wrapper
from c5_input_cleaning import input_cleaning_wrapper
from c5_input_cleaning import df

# autopep8 --in-place --aggressive --aggressive c5_pipeline.py

def c5_wrapper(dataframe):
    """
    Takes dataframe and runs it through cleaning script.
    Passes result through RF Classifier model.

    Input
    ----------
    Pandas dataframe

    Returns
    -------
    -Vector of predictions (numpy arrray)
    -Parity plot
    -Model score
    """
    # clean input dataframe
    cleaned = input_cleaning_wrapper(dataframe)

    # run through model
    classifier = rf_wrapper(cleaned)

    return classifier


print(c5_wrapper(df))
