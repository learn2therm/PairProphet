"""
This module cleans dataframe from user input for upstream classifier.
Keeps both protein sequences for reporting results.
"""

from sklearn.utils import resample
import pandas as pd

# # sample dataframe can be passed into wrapper for training
# df = pd.read_csv('learn2therm_sample_50k.csv')

# # target from Humood
# target = pd.read_csv('protein_match_50k.csv')

# # # Separate the majority and minority classes
# majority_class = target[target['hmmer_match'] is 'Yes']
# minority_class = target[target['hmmer_match'] is 'No']

# # # Undersample the majority class to match the number of minority class samples
# n_samples = len(minority_class)
# undersampled_majority = resample(
#     majority_class,
#     n_samples=n_samples,
#     replace=False)

# # # Combine the undersampled majority class with the minority class
# balanced_data = pd.concat([undersampled_majority, minority_class])

# df = pd.merge(df, target, on=['prot_pair_index'])


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
    'm_protein_len',
    't_protein_len',
    'hmmer_match',
    'structure_match',
    'norm_bit_score_m',
    'norm_bits_score_t',
    'm_protein_seq',
    't_protein_seq',
    'prot_pair_index']


def normalize_bit_scores(dataframe):
    '''
    Creates two new columns of bit score
    normalized by the protein length.

    Args:
        Dataframe (pandas dataframe)

    Returns:
        Dataframe (pandas dataframe)
    '''

    dataframe['norm_bit_score_m'] = dataframe['bit_score'] / \
        dataframe['m_protein_len']
    dataframe['norm_bit_score_t'] = dataframe['bit_score'] / \
        dataframe['t_protein_len']

    return dataframe

# need function that merges in protein_match
# need function that gets rid of Unnamed:0 and Jaccard_Score


def check_input_type(dataframe):
    '''
    Takes in input dataframe and asserts that it is the correct data type.

    Args:
        Dataframe (pandas dataframe)

    Returns:
        Dataframe (pandas dataframe)
    '''
    assert "pandas.core.frame.DataFrame" in str(
        type(dataframe)), 'Not a pandas dataframe!'

    return dataframe


def clean_input_columns(dataframe):
    '''
    Cleans out columns that are not in
    a predefined list of features.

    Args:
        Dataframe (pandas dataframe)

    Returns:
        Dataframe (pandas dataframe)
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
        Dataframe (pandas dataframe)

    Returns:
        Dataframe (pandas dataframe)
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
        Dataframe (pandas dataframe)

    Returns:
        Dataframe (pandas dataframe)
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
        Dataframe (pandas dataframe)

    Returns:
        Dataframe (pandas dataframe)
    '''
    assert 'm_protein_len' in dataframe, 'Dataframe missing mesophillic sequence!'
    assert 't_protein_len' in dataframe, 'Dataframe missing thermophillic sequence!'

    return dataframe


def input_cleaning_wrapper(dataframe):
    '''
    Takes in a pandas dataframe and runs it through each of the cleaning
    and verification steps.
    Args:
        Dataframe (pandas dataframe)

    Returns:
        Dataframe (pandas dataframe)
    '''
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
