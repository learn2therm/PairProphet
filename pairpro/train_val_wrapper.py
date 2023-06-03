"""
Wrapper functions for all of the machine learning component.
"""

from pairpro.train_val_classification import rf_wrapper
from pairpro.train_val_input_cleaning import input_cleaning_wrapper
# from pairpro.train_val_input_cleaning import df
from pairpro.train_val_featuregen import create_new_dataframe


def train_val_wrapper(dataframe, target, feature_list=None):
    '''
    Takes dataframe and runs it through cleaning script.
    Generates features with iFeatureOmegaCLI.
    Passes result through RF Classifier model.

    Args:
        Dataframe (pandas dataframe)
        Features from iFeatureOmega (list)

    Returns:
        Vector of predictions (numpy arrray)
        Parity plot
        Model score
    '''
    # clean input dataframe
    dataframe = input_cleaning_wrapper(dataframe)

    if feature_list is not None:
        # generate features from amino acid sequence
        dataframe = create_new_dataframe(dataframe, ['sequences_a.fasta', 
                                    'sequeneces_b.fasta'], 
                                    descriptors=[feature for feature in feature_list])
    else:
        pass

    #drop sequences
    dataframe.drop(columns=['m_protein_seq', 't_protein_seq'])
    # run through model
    score = rf_wrapper(dataframe, target)[1]
    model = rf_wrapper(dataframe, target)[2]

    return score, model 


