"""
The following is importable random utilites.
You will find:
- logger function
- pairwise sequence builder
"""
import pandas as pd
import itertools
import duckdb
import numpy as np
import dask.dataframe

# import distributed
import logging
logger = logging.getLogger(__name__)


def start_logger_if_necessary(logger_name: str, log_file: str, log_level,
                              filemode: str = 'a', worker: bool = False):
    """Quickly configure and return a logger that respects parallel processes.

    Parameters
    ----------
    logger_name : str
        name of logger to start or retrieve
    log_file : str
        path to file to log to
    log_level
        log level to respect
    worker: str
        name of worker using this logger
    filemode : str
        mode to apply to log file eg "a" for append
    """
    logger = logging.getLogger(logger_name)
    logger.setLevel(log_level)
    fh = logging.FileHandler(log_file, mode=filemode)
    if worker:
        # worker_name = distributed.get_worker().name
        # fh.setFormatter(logging.Formatter('%(filename)s %(worker)s - %(asctime)s %(levelname)-8s %(message)s'))
        # if len(logger.handlers) == 0:
        #     logger.addHandler(fh)
        # else:
        #     logger.handlers[-1] = fh
        # logger = logging.LoggerAdapter(logger, extra={'worker': worker_name})
        pass
    else:
        fh.setFormatter(logging.Formatter('%(filename)s - %(asctime)s %(levelname)-8s %(message)s'))
        if len(logger.handlers) == 0:
            logger.addHandler(fh)
        else:
            logger.handlers[-1] = fh
    return logger


def make_pairs(seq1_list, seq2_list, seq1_name='seq1', seq2_name='seq2',
               csv_path='./paired_seqs.csv', save=True):
    '''
    Function for building a combinatorial set of sequences from two lists.

    Args:
        seq1_list (list): List of protein sequence strings
        seq2_list (list): List of protein sequence strings
        seq1_name (str): Column name for first sequence column
        seq2_name (str): Column name for second sequence column
        csv_path (str): Path for saved .csv file
        save (bool): Saves paired sequences as .csv when True

    Returns:
        combined_df (pd.DataFrame): A dataframe with rows as all possible
        sequence pairs.
    '''
    combined = list(itertools.product(seq1_list, seq2_list))
    combined_df = pd.DataFrame(combined, columns=[seq1_name, seq2_name])

    if save is True:

        combined_df.to_csv(csv_path, index=False)

    return combined_df


# TODO: Check that skiprows fix is implemented in dask such that only first 
# partition is shortened
def split_txt(in_path, out_dir, cols, blocksize='10MB', skiprows=4, sep='\t', **kwargs):
    '''
    Converts large text files into a set of parquet files.

    Args:
        in_path (str): Path to input file.
        out_dir (str): Output directory for parquet files.
        cols (list): Desired columns from input file
        blocksize (list): Maximum size for each output file.
        skiprows (int): Number of rows to skip at beginning of text file.
        sep (str): Separator for text file.
    '''
    # Dask read_csv function is able to stream large text files to a dataframe
    # object.
    df = dask.dataframe.read_csv(f'{in_path}', names=cols, skiprows=skiprows, 
                                 blocksize=blocksize, sep=sep, header=None,
                                 **kwargs)
    df.to_parquet(out_dir)

def parqs_to_db(parq_dir, table_name, con, cols=['*']):
    '''
    Converts a DuckDB database from a collection of parquet files.

    Args:
        parq_dir (str): Path to folder containing parquet files.
        table_name (str): Name of desired table in database.
        con (duckdb.DuckDBPyConnection): Connection to existing DuckDB 
                                         database.
        cols (list): List of desired columns from parquet files.
    '''
    con.execute(f"""CREATE OR REPLACE TABLE {table_name} AS 
                    SELECT {', '.join(cols)} 
                    FROM read_parquet('{parq_dir}/*.parquet')""")
    con.commit()

def oma_sample(con, size, oma_size = 1900000000):
    '''
    Collects sample of the desired size from the OMA database.

    Args:
        con (duckdb.DuckDBPyConnection): Connection to OMA database.
        size (list): Desired number of samples.
        oma_size (int): Total size of OMA database sample is drawn from. 
                        Default 1.9 billion corresponds to total number of
                        prokaryotic pairs before cleaning.

    Returns:
        df (pd.DataFrame): A dataframe containing sample IDs and sequences.
    '''

    # Due to the size of OMA, it is necessary to sample efficiently. This is
    # done by selecting a random sample of the desired size from the 
    # prokaryotic pairs table, then joining the protein sequences to the 
    # sample.
    df = con.execute(f"""SELECT protein1_uniprot_id,
                         protein2_uniprot_id,
                         s1.sequence AS protein1_sequence,
                         s2.sequence AS protein2_sequence
                         FROM prok_pairs
                         LEFT JOIN proteins s1 
                         ON s1.id = protein1_oma_id
                         LEFT JOIN proteins s2 
                         ON s2.id = protein2_oma_id 
                         WHERE RANDOM() < {size/oma_size}""").fetchdf()
    return df