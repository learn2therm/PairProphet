"""
Wrapper functions for all of the machine learning component.
"""

from train_val_classification import rf_wrapper
from train_val_input_cleaning import input_cleaning_wrapper
from train_val_input_cleaning import df


def train_val_wrapper(dataframe):
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

print(train_val_wrapper(df))
