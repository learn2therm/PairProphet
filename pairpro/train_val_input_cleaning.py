"""
This module takes a dataframe from the data scraping component
and cleans it so that it can be passed through a machine
learning algorithm.
"""

from sklearn.utils import resample
import pandas as pd

# keep columns that can be used as features
columns_to_keep = [
    'pair_id',
    'bit_score',
    'global_gap_compressed_percent_id',
    'scaled_global_query_percent_id',
    'scaled_global_symmetric_percent_id',
    'query_align_len',
    'query_align_cov',
    'subject_align_len',
    'subject_align_cov',
    'query_len',
    'subject_len',
    'norm_bit_score_query',
    'norm_bit_score_subject',
    'query',
    'subject',
    'true_labels'
]


def verify_input_data(dataframe, blast=False, hmmer=False, structure=False):
    '''
    Verifies input dataframe and performs necessary cleaning steps.

    Args:
        dataframe: pandas DataFrame
        structure: boolean indicating whether to include 'structure_match' column

    Returns:
        dataframe: cleaned and verified pandas DataFrame
    '''

    if structure:
        columns_to_keep.append('structure_match')
    else:
        pass

    if hmmer:
        columns_to_keep.append('hmmer_match')
    else:
        pass

    # Normalize bit scores
    if blast:
        dataframe['norm_bit_score_query'] = dataframe['bit_score'] / dataframe['query_len']
        dataframe['norm_bit_score_subject'] = dataframe['bit_score'] / dataframe['subject_len']
        assert 'query_len' in dataframe.columns, 'Dataframe missing query sequence!' # Verify blast pairs
        assert 'subject_len' in dataframe.columns, 'Dataframe missing subject sequence!' # Verify blast pairs
    else:
        pass

    # Check type of dataframe
    assert isinstance(dataframe, pd.DataFrame), 'Not a pandas dataframe!'

    # Confirm data has two input sequences
    assert 'query' in dataframe.columns, 'Dataframe missing query sequence!' 
    assert 'subject' in dataframe.columns, 'Dataframe missing subject sequence!' 

    # Clean out unnecessary columns
    dataframe.drop(columns=dataframe.columns.difference(columns_to_keep), inplace=True)

    # Drop rows with NaN values
    dataframe.dropna(subset=columns_to_keep, inplace=True)

    print('The new shape of the dataframe is: {}'.format(dataframe.shape))

    return dataframe