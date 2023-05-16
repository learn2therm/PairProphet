"""
Wrapper functions for all of the machine learning component.
"""

from train_val_classification import rf_wrapper
from train_val_input_cleaning import input_cleaning_wrapper
# from train_val_input_cleaning import df
from train_val_featuregen import create_new_dataframe


def train_val_wrapper(dataframe, feature_list=None):
    """
    Takes dataframe and runs it through cleaning script.
    Generates features with iFeatureOmegaCLI.
    Passes result through RF Classifier model.

    Input
    ----------
    Pandas dataframe
    List of features from iFeatureOmega

    Returns
    -------
    -Vector of predictions (numpy arrray)
    -Parity plot
    -Model score
    """
    # clean input dataframe
    cleaned = input_cleaning_wrapper(dataframe)

    if feature_list is not None:
        # generate features from amino acid sequence
        df = create_new_dataframe(cleaned, 'meso_50k.fasta', 
                    descriptors=[feature for feature in feature_list])
    else:
        pass

    # run through model
    classifier = rf_wrapper(df)

    return classifier
