"""
To do: Raise exception for invalid inputs, try capitalization before removing
rows

OLD VERSION
# TODO: Clean up ray implementation
"""
from Bio import Align
from Bio.Align import substitution_matrices
import sys
import numpy as np
import pandas as pd
import re
import duckdb
import pickle
from joblib import Parallel, delayed
from joblib.externals.loky import get_reusable_executor




class PicklablePairwiseAligner(Align.PairwiseAligner):
    """
    Picklable version of the PairwiseAligner class from the Bio.Align module.
    """
    def __getstate__(self):
        return
    def __setstate__(self, state):
        return

def make_blast_df(df_in, cpus=2, mode='local', path='./data/blast_db.db', module_path='.'):
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

    Notes:
        This function uses joblib to parallelize the alignment process. The
        number of cores used is determined by the cpus parameter. The 
        alignment scores are calculated using the BLOSUM62 substitution matrix.
    """
    
    # Rename input data columns for compatibility
    df = df_in.rename(columns={'protein1_sequence': 'query', 'protein2_sequence': 'subject'})


    # Remove any rows with NaN or containing non-amino acid letters
    df = df.dropna().reset_index(drop=True)

    # All valid 1-letter amino acid codes.
    amino_acids = set("CSTAGPDEQNHRKMILVWYF")

    invalid_rows = []
    for index, row in df.iterrows():
        if not set(row['query']).issubset(amino_acids) or not set(row['subject']).issubset(amino_acids):
            invalid_rows.append(index)
    
    if invalid_rows:
        print(f"Found and skipped {len(invalid_rows)} invalid row(s) containing invalid amino acid sequences.")
        df = df.drop(index=invalid_rows).reset_index(drop=True)

    # Metrics to be calculated
    metrics = ['local_gap_compressed_percent_id',
               'scaled_local_query_percent_id',
               'scaled_local_symmetric_percent_id',
               'query_align_len',
               'query_align_cov',
               'subject_align_len',
               'subject_align_cov',
               'bit_score']

    # Initialize PairWiseAligner with required parameters for model
    aligner = PicklablePairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = mode

    # Iterate through all pairs and calculate best alignment and metrics
    cols = metrics + ['pair_id', 'protein1_uniprot_id', 'protein2_uniprot_id']
    final_data = pd.DataFrame(columns=cols)

    def aligner_align(row, aligner):

        # Initialize PairWiseAligner with required parameters for model locally (messy code)
        aligner = PicklablePairwiseAligner()
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score = -11
        aligner.extend_gap_score = -1
        aligner.mode = mode
  
        index, data = row  # Unpack the tuple into index and data (Series)
        subject = data['subject']
        query = data['query']

        alignment = aligner.align(subject, query)
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

        # Collect calculated metrics into a list for addition to final_data
        new_row = [gap_comp_pct_id, scaled_local_query_percent_id,
                scaled_local_symmetric_percent_id, query_length,
                query_cov, subject_length, subject_cov, best_alignment.score,
                data['pair_id'], data['protein1_uniprot_id'], data['protein2_uniprot_id']]
        return new_row
    
    final_data = Parallel(n_jobs=cpus)(delayed(aligner_align)(row, aligner) for row in df.iterrows())

    # Convert the list of results to a DataFrame
    blast_df = pd.DataFrame(final_data, columns=cols)

    # Merge final_df with sequences and ids from input df
    # blast_df = df.merge(final_data_df, on=['pair_id', 'protein1_uniprot_id', 'protein2_uniprot_id'])
    con = duckdb.connect(path)
    cmd = """CREATE OR REPLACE TABLE protein_pairs
             AS SELECT * FROM blast_df"""
    con.execute(cmd)

    return blast_df, con

def blast_pairs(df_in, cpus=2, mode='local'):
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

    Notes:
        This function uses joblib to parallelize the alignment process. The
        number of cores used is determined by the cpus parameter. The 
        alignment scores are calculated using the BLOSUM62 substitution matrix.
    """
    get_reusable_executor().shutdown(wait=True)
    # Remove any rows with NaN or containing non-amino acid letters
    df = df_in.dropna().reset_index(drop=True)

    # All valid 1-letter amino acid codes.
    amino_acids = set("CSTAGPDEQNHRKMILVWYF")

    invalid_rows = []
    for index, row in df.iterrows():
        if not set(row['query']).issubset(amino_acids) or not set(row['subject']).issubset(amino_acids):
            invalid_rows.append(index)
    
    if invalid_rows:
        print(f"Found and skipped {len(invalid_rows)} invalid row(s) containing invalid amino acid sequences.")
        df = df.drop(index=invalid_rows).reset_index(drop=True)

    # Metrics to be calculated
    metrics = ['local_gap_compressed_percent_id',
               'scaled_local_query_percent_id',
               'scaled_local_symmetric_percent_id',
               'query_align_len',
               'query_align_cov',
               'subject_align_len',
               'subject_align_cov',
               'bit_score']

    # Initialize PairWiseAligner with required parameters for model
    aligner = PicklablePairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = mode

    # Iterate through all pairs and calculate best alignment and metrics
    cols = metrics + ['pair_id', 'query_id', 'subject_id']
    final_data = pd.DataFrame(columns=cols)

    def aligner_align(row):

        try:
            # Initialize PairWiseAligner with required parameters for model locally (messy code)
            aligner = PicklablePairwiseAligner()
            aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
            aligner.open_gap_score = -11
            aligner.extend_gap_score = -1
            aligner.mode = mode
    
            # index, data = row  # Unpack the tuple into index and data (Series)
            subject = row['subject']
            query = row['query']

            alignment = aligner.align(subject, query)
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

            # Collect calculated metrics into a list for addition to final_data
            new_row = [gap_comp_pct_id, scaled_local_query_percent_id,
                    scaled_local_symmetric_percent_id, query_length,
                    query_cov, subject_length, subject_cov, best_alignment.score,
                    row['pair_id'], row['query_id'], row['subject_id']]
            return new_row
        except Exception as e:
            print(f"Error in row: {row}: {e}")
            return [None]*len(cols)
    
    final_data = Parallel(n_jobs=cpus, verbose=50, timeout=2000)(delayed(aligner_align)((index, row)) for _, row in df.iterrows())

    # Convert the list of results to a DataFrame
    blast_df = pd.DataFrame(final_data, columns=cols)

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
