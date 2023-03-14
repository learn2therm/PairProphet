"""
This module takes a dataframe from the data scraping component
and cleans it so that it can be passed through a machine
learning algorithm.
"""

# import pandas as pd

# sample dataframe can be passed into wrapper. Commented out for now
# df = pd.read_csv('learn2therm_sample_50k.csv')

# df['protein_match'] = df['t_protein_desc'] == df['m_protein_desc']

# keep columns we are interested in
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


def check_input_type(dataframe):
    """
    Takes in input dataframe and asserts that it is the correct data type.
    """
    assert "pandas.core.frame.DataFrame" in str(
        type(dataframe)), 'Not a pandas dataframe!'

    return dataframe


def clean_input_columns(dataframe):
    """
    We want to clean certain columns out of the Pfam dataframe.
    Need to eliminate identifier columns + columns that don't have
    relationship with the target.

    Input: Pandas dataframe (from Pfam)
    Output: Updated dataframe.
    """

    for title in dataframe:
        if title not in columns_to_keep:
            dataframe = dataframe.drop(columns=title)
        else:
            pass

    return dataframe


def verify_input_columns(dataframe):
    """
    This function raises an error is one of the columns we need for the model is not
    present in the dataframe.

    Input: Pandas dataframe.
    Output: Pandas dataframe.
    """
    for title in columns_to_keep:

        if title not in dataframe:
            raise KeyError
        else:
            pass

    return dataframe


def check_input_nans(dataframe):
    """
    Checks for NaN values in input dataframe. Removes rows with NaN values present.

    Input: Pandas dataframe
    Output: Pandas dataframe

    """
    has_nan = dataframe.isna().any().any()
    nan_rows = dataframe[dataframe.isna().any(axis=1)]

    if has_nan:
        print('Dataframe has {} rows with NaN values!'.format(len(nan_rows)))
    else:
        print("DataFrame does not have any NaN values.")

    # Drop rows with NaN's
    dataframe = dataframe.dropna()
    print('Dataframe now has {} rows.'.format(len(dataframe)))

    return dataframe


def verify_protein_pairs(dataframe):
    """
    Checks that input data has two protein sequences. Will need to generalize this function other data sets
    to simply make sure two sequences are entered. Code below is for our protein database
    """
    assert 'm_protein_len' in dataframe, 'Dataframe missing mesophillic sequence!'
    assert 't_protein_len' in dataframe, 'Dataframe missing thermophillic sequence!'

    print('OK!')
    return dataframe


def input_cleaning_wrapper(dataframe):
    """
    Takes in a pandas dataframe and runs it through each of the cleaning
    and verification steps.

    Input: Pandas dataframe
    Output: Pandas dataframe

    """

    # check type of dataframe
    check = check_input_type(dataframe)

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
