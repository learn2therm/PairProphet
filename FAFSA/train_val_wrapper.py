"""
Wrapper functions for all of the machine learning component.
"""

from train_val_classification import rf_wrapper
from train_val_input_cleaning import input_cleaning_wrapper
from train_val_input_cleaning import df
from train_val_featuregen import create_new_dataframe


def train_val_wrapper(dataframe):
    """
    Takes dataframe and runs it through cleaning script.
    Generates features with iFeatureOmegaCLI.
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

    #generate features from amino acid sequence
    df_with_features = create_new_dataframe(cleaned, 'meso_50k.fasta', descriptors=['AAC', 'GAAC'])

    # run through model
    classifier = rf_wrapper(df_with_features)

    return classifier

print(train_val_wrapper(df))
