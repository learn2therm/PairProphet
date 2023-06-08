"""
This module cleans dataframe from user input for upstream classifier.
Keeps both protein sequences for reporting results.
"""

from sklearn.utils import resample
import pandas as pd


# keep columns that can be used as features
columns_to_keep = [
    'bit_score',
    'local_gap_compressed_percent_id',
    'scaled_local_query_percent_id',
    'scaled_local_symmetric_percent_id',
    'query_align_len',
    'query_align_cov',
    'subject_align_len',
    'subject_align_cov',
    'query_len',
    'subject_len',
    'hmmer_match',
    'norm_bit_score_query',
    'norm_bit_score_subject',
    'query',
    'subject'
]


def normalize_bit_scores(dataframe):
    '''
    Creates two new columns of bit score
    normalized by the protein length.

    Args:
        pandas dataframe

    Returns:
        pandas dataframe
    '''
    dataframe['norm_bit_score_query'] = dataframe['bit_score'] / \
        dataframe['query_len']
    dataframe['norm_bit_score_subject'] = dataframe['bit_score'] / \
        dataframe['subject_len']

    return dataframe


def check_input_type(dataframe):
    '''
    Takes in input dataframe and asserts that it is the correct data type.

    Args:
        pandas dataframe

    Returns:
        pandas dataframe
    '''
    assert "pandas.core.frame.DataFrame" in str(
        type(dataframe)), 'Not a pandas dataframe!'

    return dataframe


def clean_input_columns(dataframe):
    '''
    Cleans out columns that are not in
    a predefined list of features.

    Args:
        pandas dataframe

    Returns:
        pandas dataframe
    '''
    for title in dataframe:
        if title not in columns_to_keep:
            dataframe = dataframe.drop(columns=title)
        else:
            pass

    return dataframe


def verify_input_columns(dataframe):
    '''
    Asserts that columns we want to keep
    remain in the dataframe.

    Args:
        pandas dataframe

    Returns:
        pandas dataframe
    '''
    for title in columns_to_keep:

        if title not in dataframe:
            raise KeyError
        else:
            pass

    return dataframe


def check_input_nans(dataframe):
    '''
    Checks for NaN values in input dataframe.
    Removes rows with NaN values present.
    Args:
        pandas dataframe

    Returns:
        pandas dataframe
    '''
    has_nan = dataframe.isna().any().any()
    nan_rows = dataframe[dataframe.isna().any(axis=1)]

    if has_nan:
        print('Dataframe has {} rows with NaN values!'.format(len(nan_rows)))
    else:
        print("DataFrame does not have any NaN values.")

    # Drop rows with NaN's
    dataframe = dataframe.dropna()

    return dataframe


def verify_protein_pairs(dataframe):
    '''
    Checks that input data has two protein sequences
    with simple assert statements.
    Args:
        pandas dataframe

    Returns:
        pandas dataframe
    '''
    assert 'query_len' in dataframe, 'Dataframe missing query sequence!'
    assert 'subject_len' in dataframe, 'Dataframe missing subject sequence!'

    return dataframe


def input_cleaning_wrapper(dataframe, structure):
    '''
    Takes in a pandas dataframe and runs it through each of the cleaning
    and verification steps.
    Args:
        pandas dataframe

    Returns:
        pandas dataframe
    '''
    if structure:
        columns_to_keep.append('structure_match')
    else:
        pass

    # normalize bit scores
    normed = normalize_bit_scores(dataframe)

    # check type of dataframe
    check = check_input_type(normed)

    # clean out unnecessary columns
    clean = clean_input_columns(check)

    # verify necessary columns are present
    verify_input = verify_input_columns(clean)

    # check for NaN's
    check_nans = check_input_nans(verify_input)

    # verify every protein has a pair
    verify_pairs = verify_protein_pairs(check_nans)

    print('The new shape of the dataframe is:{}'.format(verify_pairs.shape))

    return verify_pairs
