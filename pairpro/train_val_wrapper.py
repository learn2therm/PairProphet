"""
Wrapper functions for all of the machine learning component.
"""

from pairpro.train_val_classification import rf_wrapper
from pairpro.train_val_input_cleaning import verify_input_data
from pairpro.train_val_featuregen import create_new_dataframe


def train_val_wrapper(dataframe, target, blast=False, 
                      hmmer=False, structure=False, features=False):
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
    dataframe = verify_input_data(dataframe, blast, hmmer, structure)

    if features is True:
        feature_list = [
            'AAC', 'GAAC', 'DistancePair', 'CTDC', 'CTDT', 'CTDD', 'CTriad', 'GDPC type 1', 'GDPC type 2', 'CKSAAGP type 1', 'CKSAAGP type 2',
            'PseKRAAC type 2', 'PseKRAAC type 3A', 'PseKRAAC type 7', 'PseKRAAC type 9', 'Geary', 'APAAC', 'QSOrder']
        # generate features from amino acid sequence
        dataframe = create_new_dataframe(dataframe, ['sequences_a.fasta',
                                                     'sequeneces_b.fasta'],
                                         descriptors=[feature for feature in feature_list])
    else:
        pass

    # drop sequences
    dataframe.drop(
        columns=[
            'subject',
            'query',
            'pair_id'],
        inplace=True)
    # run through model
    score, model = rf_wrapper(dataframe, target)

    return score, model
