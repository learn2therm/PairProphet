"""
The script performs search on Protein Data Bank (PDB) and return output.
The goal is to return proteins with structural similarity to query.
Current code allows users to access PDB, but only returns proteins with sequence homology and not structural similarity.
We plan to improve this in Spring 2023.

The packages you need to run this script are the following:
- pyPDB
- pandas
"""

import pandas as pd
import pypdb
import json

def scrapepdb(seq: str):
    """
    ------------
    Parameters
    ------------
    seq: str
        a str that represents the amino acid sequence.
    
    Returns:
    ------------
    df_result_reind : pandas.core.frame.DataFrame
        a dataframe containing identified proteins with sequence similarity to query. Percent identity, bitscore and e-value are included for evaluation.
    """

    # Search the PDB using the provided sequence. 
    return_polymer_entity = pypdb.Query(seq, query_type='sequence', return_type='polymer_entity')
    # Drop unnecessary information (result_type: polymer_entity, total_count) to format it into DataFrame.
    result_set = return_polymer_entity.search()['result_set']
    # Sort the result_set into DataFrame.
    df_result = pd.json_normalize(result_set, ['services', 'nodes', 'match_context'], meta=['identifier', 'score', ['services', 'service_type']])
    # Order and keep certain columns based on what we care most, also this can be modified anytime.
    df_result_reind = df_result.reindex(columns=['identifier', 'score', 'sequence_identity', 'evalue', 'bitscore'])
    return df_result_reind