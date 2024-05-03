"""
To do: Raise exception for invalid inputs, try capitalization before removing
rows

OLD VERSION
# TODO: Clean up ray implementation
"""
from Bio import Align
from Bio.Align import substitution_matrices
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
import functools
import sys
import multiprocessing
import numpy as np
import pandas as pd
import re
import duckdb
import pickle
from joblib import Parallel, delayed
from joblib.externals.loky import get_reusable_executor


# Setup basic configuration for logging
logger = logging.getLogger(__name__)
# logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')



class PicklablePairwiseAligner(Align.PairwiseAligner):
    """
    Picklable version of the PairwiseAligner class from the Bio.Align module.
    """
    def __getstate__(self):
        return
    def __setstate__(self, state):
        return

# Global variable in the worker process
global_aligner = None
    
def worker_init():
    """
    Initializes the aligner with the BLOSUM62 substitution matrix.
    """
    global global_aligner
    global_aligner = PicklablePairwiseAligner()
    global_aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    global_aligner.open_gap_score = -11
    global_aligner.extend_gap_score = -1
    global_aligner.mode = 'global' # 'local' or 'global'


def process_row_nuprok(row):
    """
    Processes a single row of a DataFrame to calculate alignment metrics.
    Nuprok version.
    """
    try:
        subject = row['subject']
        query = row['query']

        alignment = global_aligner.align(subject, query)
        best_alignment = max(alignment, key=lambda x: x.score)

        alignment_str = format(best_alignment)
        alignment_lines = alignment_str.split('\n')

        seq1_aligned = alignment_lines[0]
        seq2_aligned = alignment_lines[2]

        # Coverage and sequence length
        subject_cov = sum(c != '-' for c in seq1_aligned)
        query_cov = sum(c != '-' for c in seq2_aligned)
        subject_length = len(seq1_aligned)
        query_length = len(seq2_aligned)
        subject_cov /= subject_length
        query_cov /= query_length

        # Percent ids
        n_matches, n_gaps, n_columns, n_comp_gaps = get_matches_gaps(
            seq2_aligned, seq1_aligned)
        gap_comp_pct_id = gap_compressed_percent_id(
            n_matches, n_gaps, n_columns, n_comp_gaps)
        scaled_local_symmetric_percent_id = 2*n_matches / (subject_length + query_length)
        scaled_local_query_percent_id = n_matches / query_length

        # Collect calculated metrics
        new_row = {
            'pair_id': row.get('pair_id'),
            'protein1_uniprot_id': row.get('protein1_uniprot_id'),
            'protein2_uniprot_id': row.get('protein2_uniprot_id'),
            'bit_score': best_alignment.score,
            'local_gap_compressed_percent_id': gap_comp_pct_id,
            'scaled_local_query_percent_id': scaled_local_query_percent_id,
            'scaled_local_symmetric_percent_id': scaled_local_symmetric_percent_id,
            'query_align_len': query_length,
            'query_align_cov': query_cov,
            'subject_align_len': subject_length,
            'subject_align_cov': subject_cov,
            'query_len': len(query),
            'subject_len': len(subject),
        }
        
        return new_row
    except Exception as e:
        return {'pair_id': row.get('pair_id'), 'details': str(e)}
    

def process_row(row):
    """
    Processes a single row of a DataFrame to calculate alignment metrics.
    General version.
    """
    try:
        subject = row['subject']
        query = row['query']

        alignment = global_aligner.align(subject, query)
        best_alignment = max(alignment, key=lambda x: x.score)

        alignment_str = format(best_alignment)
        alignment_lines = alignment_str.split('\n')

        seq1_aligned = alignment_lines[0]
        seq2_aligned = alignment_lines[2]

        # Coverage and sequence length
        subject_cov = sum(c != '-' for c in seq1_aligned)
        query_cov = sum(c != '-' for c in seq2_aligned)
        subject_length = len(seq1_aligned)
        query_length = len(seq2_aligned)
        subject_cov /= subject_length
        query_cov /= query_length

        # Percent ids
        n_matches, n_gaps, n_columns, n_comp_gaps = get_matches_gaps(
            seq2_aligned, seq1_aligned)
        gap_comp_pct_id = gap_compressed_percent_id(
            n_matches, n_gaps, n_columns, n_comp_gaps)
        scaled_local_symmetric_percent_id = 2*n_matches / (subject_length + query_length)
        scaled_local_query_percent_id = n_matches / query_length

        # Collect calculated metrics
        new_row = {
            'pair_id': row.get('pair_id'),
            'query_id': row.get('query_id'),
            'subject_id': row.get('subject_id'),
            'bit_score': best_alignment.score,
            'local_gap_compressed_percent_id': gap_comp_pct_id,
            'scaled_local_query_percent_id': scaled_local_query_percent_id,
            'scaled_local_symmetric_percent_id': scaled_local_symmetric_percent_id,
            'query_align_len': query_length,
            'query_align_cov': query_cov,
            'subject_align_len': subject_length,
            'subject_align_cov': subject_cov,
            'query_len': len(query),
            'subject_len': len(subject),
        }
        
        return new_row
    except Exception as e:
        return {'pair_id': row.get('pair_id'), 'details': str(e)}

def preprocess_alignment_dataframe(df_in):
    """
    Preprocesses the input dataframe by removing rows with NaN values and
    sequences containing non-amino acid letters.

    Args:
        df_in (pandas.core.DataFrame): A dataframe that should contain protein pair sequences
                                       and their associated ids.
    
    Returns:
        df (pandas.core.DataFrame): A cleaned dataframe with the input sequences and ids.
    """
    # Remove any rows with NaN or containing non-amino acid letters
    df = df_in.dropna().reset_index(drop=True)

    # All valid 1-letter amino acid codes.
    amino_acids = set("CSTAGPDEQNHRKMILVWYF")

    # Check if all sequences are valid
    valid_query = df['query'].apply(lambda x: set(x).issubset(amino_acids))
    valid_subject = df['subject'].apply(lambda x: set(x).issubset(amino_acids))
    valid_rows = valid_query & valid_subject

    # Filter out invalid rows
    df = df[valid_rows].reset_index(drop=True)
    print(f"Found and skipped {len(df_in) - len(df)} invalid row(s) containing invalid amino acid sequences.")

    return df

def alignment_worker_og(chunk):
    """
    Processes a chunk of protein sequence pairs to calculate alignment metrics,
    including returning identifiers from the original DataFrame.
    """
    logging.debug("Starting alignment for a new chunk with %d rows,", len(chunk))
    # Apply process_row_nuprok to each row in the chunk
    results = chunk.apply(process_row_nuprok, axis=1, result_type='expand')
    return pd.DataFrame(results)    


def make_blast_df(df_in, cpus=2, path='./data/blast_db.db'):
    """
    This function generates pairwise alignment scores for a set of protein
    sequences.

    Args:
        df (pandas.core.DataFrame): A dataframe that should contain protein pair sequences
                                    and their associated ids.
        cpus (int): The number of cores to use for parallelization. Default: 2.
        mode (str): Alignment type is 'local' or 'global'. Default: 'local'.

    Returns:
        blast_df (pandas.core.DataFrame): A dataframe with the input sequence
                                          pairs, associated id values, and
                                          alignment scores.
        {optional} con (duckdb.DuckDB): A connection to the DuckDB database. 

    Notes:
        This function parallelizes the alignment process across multiple cores
        using the multiprocessing module. The number of cores (cpus) will always
        be one less than the total number of available cores to prevent
        freezing the system.
    """
    logging.debug("Starting make_blast_df with %d CPUs", cpus)
    # Rename input data columns for compatibility
    df = df_in.rename(columns={'protein1_sequence': 'query', 'protein2_sequence': 'subject'})


    # Remove any rows with NaN or containing non-amino acid letters
    df = df.dropna().reset_index(drop=True)

    # All valid 1-letter amino acid codes.
    amino_acids = set("CSTAGPDEQNHRKMILVWYF")

    df = df[df['query'].apply(lambda x: set(x).issubset(amino_acids)) & df['subject'].apply(lambda x: set(x).issubset(amino_acids))]

    if len(df) != len(df_in):
        print(f"Filtered out {len(df_in) - len(df)} rows with invalid amino acid sequences.")
    
    # Split the dataframe into chunks
    chunks = np.array_split(df, cpus-1)
    chunk_tuples = [(chunk,) for chunk in chunks]  # Wrap each chunk in a tuple

    # Pool to procces each chunk in parallel
    with multiprocessing.Pool(processes=cpus-1, initializer=worker_init) as pool:
        logging.debug("Multiprocessing pool opened with %d processes", cpus)
        results = pool.starmap(alignment_worker_og, chunk_tuples)
        logging.debug("Data processing complete, closing multiprocessing pool")

    # Concatenate the results into a single dataframe
    blast_df = pd.concat(results, ignore_index=True)

    logging.debug("make_blast_df completed, returning DataFrame")

    # Store the resulting DataFrame in a DuckDB database
    store_in_duckdb(blast_df, path)

    return blast_df

def store_in_duckdb(blast_df, path):
    """
    Stores the resulting DataFrame in a DuckDB database.
    
    Args:
        blast_df (pandas.core.DataFrame): A DataFrame containing the pairwise
                                          alignment results.
        path (str): The path to the DuckDB database file.
    """
    # create a duckdb database and store the blast_df
    con = duckdb.connect(path)
    cmd = """CREATE OR REPLACE TABLE protein_pairs
             AS SELECT * FROM blast_df"""
    con.execute(cmd)
    con.close()


def make_blast_df_sequential(df_in, path='./data/blast_db.db'):
    """
    Test function to run the alignment sequentially.
    """
    df = df_in.rename(columns={'protein1_sequence': 'query', 'protein2_sequence': 'subject'})
    df = df.dropna().reset_index(drop=True)
    amino_acids = set("CSTAGPDEQNHRKMILVWYF")
    df = df[df['query'].apply(lambda x: set(x).issubset(amino_acids)) & df['subject'].apply(lambda x: set(x).issubset(amino_acids))]
    
    aligner_params = {
        'substitution_matrix': substitution_matrices.load("BLOSUM62"),
        'open_gap_score': -11,
        'extend_gap_score': -1,
        'mode': 'global'
    }

    results = []
    for chunk in np.array_split(df, 1):  # Essentially not splitting the DataFrame
        result = alignment_worker_og(chunk, aligner_params)
        results.append(result)

    blast_df = pd.concat(results, ignore_index=True)
    store_in_duckdb(blast_df, path)
    return blast_df


def alignment_worker(chunk):
    """
    Processes a chunk of protein sequence pairs to calculate alignment metrics,
    including returning identifiers from the original DataFrame.
    """
    logging.debug("Starting alignment for a new chunk with %d rows,", len(chunk))
    # Apply process_row to each row in the chunk
    results = chunk.apply(process_row, axis=1, result_type='expand')
    return pd.DataFrame(results)    


def blast_pairs(df_in, cpus=2):
    """
    This function generates pairwise alignment scores for a set of protein
    sequences.

    Args:
        df (pandas.core.DataFrame): A dataframe that should contain protein pair sequences
                                    and their associated ids.
        cpus (int): The number of cores to use for parallelization. Default: 2.

    Returns:
        blast_df (pandas.core.DataFrame): A dataframe with the input sequence
                                          pairs, associated id values, and
                                          alignment scores.

    Notes:
        This function parallelizes the alignment process across multiple cores
        using the multiprocessing module. The number of cores (cpus) will always
        be one less than the total number of available cores to prevent
        freezing the system.
    """
    # Preprocess the input dataframe
    df_cleaned = preprocess_alignment_dataframe(df_in)

    # Split the dataframe into chunks
    chunks = np.array_split(df_cleaned, cpus-1)
    chunk_tuples = [(chunk,) for chunk in chunks]  # Wrap each chunk in a tuple

    # Pool to procces each chunk in parallel
    with multiprocessing.Pool(processes=cpus-1, initializer=worker_init) as pool:
        logging.debug("Multiprocessing pool opened with %d processes", cpus)
        results = pool.starmap(alignment_worker_og, chunk_tuples)
        logging.debug("Data processing complete, closing multiprocessing pool")

    # Concatenate the results into a single dataframe
    blast_df = pd.concat(results, ignore_index=True)

    logging.debug("make_blast_df completed, returning DataFrame")

    return blast_df


def get_matches_gaps(query, subject):
    """
    Parses sequence alignment text to calculate the number of matches, gaps,
    compressed gaps, and total columns.

    Args:
        query (str): Query aligned sequence.
        subject (str): Subject aligned sequence.

    Returns:
        n_matches (int): Number of matching amino acids in the sequence
                         alignment.
        n_gaps (int): Total number of gaps across both aligned sequences.
        n_columns (int): Length of the aligned query sequence.
        n_comp_gaps (int): Number of compressed gaps.
    """
    n_columns = len(query)
    n_gaps = sum(i == '-' for i in query) + sum(i == '-' for i in subject)
    n_matches = sum(i == j and i != '-' for i, j in zip(query, subject))
    pattern = r'(-{1,})'

    instances = re.findall(pattern, query) + re.findall(pattern, subject)
    n_comp_gaps = len(instances)

    return n_matches, n_gaps, n_columns, n_comp_gaps


def gap_compressed_percent_id(n_matches, n_gaps, n_columns, n_comp_gaps):
    """
    Calculates the percent id with compressed gaps.

    Args:
        n_matches (int): Number of matches in match columns
        n_gaps (int): Number of gaps in match columns
        n_columns (int): Total number of alignment match columns
        n_compressed_gaps (int): Number of compressed gaps in match columns

    Returns:
        n_matches / (n_columns - n_gaps + n_comp_gaps)
    """
    return n_matches / (n_columns - n_gaps + n_comp_gaps)


def sequence_validate(seq, alph):
    """
    Makes sure sequence complies with alphabet.

    Args:
        seq (int): Number of matches in match columns
        alph (int): Number of gaps in match columns

    Returns:
        (bool): True if sequence is valid, False if not
    """
    extra = set(seq) - alph
    return not extra
