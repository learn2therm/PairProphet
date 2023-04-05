"""
This script takes a user defined dataframe and an integer k, which send HTTPs requests to the HMMER API
the packages you need to run this script are:
- pandas
- requests
- urllib.parse
- time
"""

import pandas as pd
import requests
import urllib.parse
import time


def hmmerscanner(df: pd.DataFrame, k: int):
    """
    This function sends HTTP requests to the HMMER API to get information for protein sequences.
    -------------
    Parameters:
    -------------
    df: pandas.core.DataFrame
        A DataFrame that has string amino acid sequences. This function has been used Meso s
        equence, we can change that to Thermo sequence according to our needed.
    k: int
        The number of sequences to scan.
    -------------
    Raises:
    -------------
    Exception:
        Raises an exception if the status is pending for too long, if the internet isn't working,
        or if the URL system doesn't wholly answer.
    -------------
    Returns:
    -------------
    results_df: pandas.core.DataFrame
        All the families are in the rows, and we have many columns that show the information that
        we need in the future. We can drop some columns and keep the needed information.
    """
    # Check if we need to use the local function instead of the API for large values of k.
    if k > 300:
        print("Use local function for the number of sequences more than 300.")
        return pd.DataFrame()

    # Create an empty DataFrame to store the results.
    results_df = pd.DataFrame()

    # Loop through the sequences to check them.
    for i in range(k):
        # This is for meso protein sequences; we can change that in the future according to our request.
        sequence = df['m_protein_seq'][i]

        # Send an HTTP request to the HMMER API to get information for the current sequence.
        url = 'https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan'
        headers = {'Content-Type': 'application/x-www-form-urlencoded',
                   'Accept': 'application/json'}
        data = {'hmmdb': 'pfam', 'seq': f'>seq\n{sequence}'}
        data = urllib.parse.urlencode(data).encode('ascii')
        response = requests.post(url, headers=headers,
                                 data=data, allow_redirects=False)
        redirect_url = response.headers.get('Location')

        if redirect_url is None:
            # If the server doesn't work, show this error.
            print("Error: No redirect URL found in response.")
        elif redirect_url == 'late':
            # Raises an exception if the status is pending for too long.
            response.raise_for_status()
            time.sleep(180)
            raise IOError("Error notice after 3 minutes.")
        else:
            response2 = requests.get(redirect_url, headers=headers)

            # Put the results in the empty DataFrame.
            results = response2.json()
            hits = results['results']['hits']
            dfff = pd.json_normalize(
                hits, 'domains', ['acc', 'name', 'score', 'evalue', 'pvalue', 'desc'])
            dfff.insert(0, 'sequence', sequence)
            dfff = dfff.set_index('sequence')
            results_df = pd.concat([results_df, dfff])
            if redirect_url == 'late':
                # Raises an exception if the status is pending for too long.
                response2.raise_for_status()
                time.sleep(180)
                raise IOError("Error notice after 3 minutes.")

    return results_df