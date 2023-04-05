'''
This package imports data to be used in the FAFSA model. Current support for DuckDB objects
generated from upstream component and for user-generated DataFrames with correct structure.
Future updates will add option to calculate alignment features during dataset assembly, but
current version requires user to supply them.

Functions:
    fetch_data: Generates input data for the FAFSA model in the form of a DataFrame with
    alignment features as columns. User can control size of sample and sampling method.
'''
import os
import sys

import pandas as pd

module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)
    
from FAFSA.preprocessing import connect_db

def fetch_data(path, form = 'duckdb', size: int = 1000, method = 'random', idx_range: list = [0,0],
               chunksize: int = 10 ** 5):
    '''
    Pulls data from DuckDB database or pandas DataFrame for input to FAFSA model.

    Args:
        path (str): Path to database or DataFrame input.
        form (str): Identifies import method for either DuckDB database or DataFrame.
        size (int): Sample size to be passed to FAFSA model. Superceded by idx_range when using sequential sampling.
        method (str): Specifies random, chunked, or numeric sampling.
        idx_range (list): Only applies to numeric sampling. Min and max index range for numeric
                          sampling.
        chunksize (int): Only applies to chunk sampling for large datasets. Sets size of each chunk.

    Returns:
        fafsa_df (pandas.DataFrame): DataFrame formatted for FAFSA model.

    Raises:
        ValueError: form must be one of ['csv', 'duckdb']. Other data types not yet supported.
        ValueError: method must be one of ['random', 'numeric', 'chunk']
        ValueError: Index range must be positive.
        ValueError: Index range must be contained in input data.
    '''

    forms = ['csv', 'duckdb']
    methods = ['random', 'numeric', 'chunk']

    if form not in forms:
        raise ValueError(f'Invalid argument passed to form. Expected one of: {forms}')

    if method not in methods:
        raise ValueError(f'Invalid argument passed to method. Expected one of: {methods}')

    if idx_range[0] > idx_range[1]:
        raise ValueError('Index range must be positive.')
        
    if idx_range[0] < 0:
        raise ValueError('Indices must be positive.')
        
    if form == 'csv':

        # Uses pandas random sampling of DataFrame.
        if method == 'random':

            fafsa_df = pd.read_csv(path).sample(size)

        # Uses pandas chunking of DataFrame with final concatention.
        elif method == 'chunk':

            dfs = pd.read_csv(path, nrows = size, chunksize=chunksize)
            fafsa_df = pd.concat(dfs, ignore_index = True)

        # Selects only rows specified by user.
        else:

            if idx_range == [0,0]:
                idx_range = [0, size]
           
            try:
                fafsa_df = pd.read_csv(path, skiprows = lambda x: x not in range(idx_range[0],
                                                                                 idx_range[1]+1))

            except:
                raise ValueError('idx_range not in input data.')

    if form == 'duckdb':

        con = connect_db(path)

        # Uses duckdb random sampling via SQL.
        if method == 'random':

            sample_cmd = f"""SELECT *
                             FROM fafsa_final
                             USING SAMPLE {size}"""
            fafsa_df = con.execute(sample_cmd).df()

        # Performs a series of SQL queries until size is met.
        elif method == 'chunk':

            offset = 0
            dfs = []

            while True:
                
                # Checks if size has been exceeded. Sets to size if it has.
                if offset + chunksize > size:
                    chunksize = size - offset

                chunk_cmd = f"""SELECT *
                                FROM fafsa_final
                                LIMIT {chunksize}
                                OFFSET {offset}
                                ORDER BY prot_pair_index"""
                dfs.append(con.execute(chunk_cmd).df())
                offset += chunksize

                if offset > size:
                    break

                fafsa_df = pd.concat(dfs, ignore_index = True)

        # Selects only rows specified by the user.
        else:

            if idx_range == [0,0]:
                idx_range = [0, size]

            num_cmd = f"""SELECT *
                         FROM fafsa_final
                         LIMIT {idx_range[1] - idx_range[0]}
                         OFFSET {idx_range[0]}
                         ORDER BY prot_pair_index"""
            try:
                fafsa_df = con.execute(num_cmd).df()
                
            except:
                raise ValueError('idx_range not in input data.')
                
    return fafsa_df
