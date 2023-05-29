from Bio import Align
from Bio.Align import substitution_matrices
import re
import numpy as np
import pandas as pd


def make_blast_df(df, mode='local'):
    """
    This function generates pairwise alignment scores for a set of protein
    sequences.

    Args:
        df (pandas.core.DataFrame): A 2-column DataFrame containing the query
                                    and subject sequences for alignment.
        mode (str): Alignment type is 'local' or 'global'. Default: 'local'.

    Returns:
        blast_df (pandas.core.DataFrame): A dataframe with the input sequence
                                          pairs, associated id values, and
                                          alignment scores.
    """
    # Rename input data columns for compatibility
    original_cols = df.columns
    df.rename(columns={original_cols[0]: 'query', original_cols[1]: 'subject'},
              inplace=True)

    # Generate protein ids for each unique query and subject sequence
    n_unique_1 = np.unique(df['query'])
    n_unique_2 = np.unique(df['subject'])

    qid_dict = dict(zip(n_unique_1, range(len(n_unique_1))))
    sid_dict = dict(zip(n_unique_2, range(len(n_unique_2))))

    df['qid'] = [qid_dict[i] for i in df.iloc[:, 0]]
    df['sid'] = [sid_dict[i] for i in df.iloc[:, 1]]

    # Use unique row index as a pair id
    df['pid'] = df.index

    # Metrics to be calculated
    metrics = ['local_gap_compressed_percent_id',
               'scaled_local_query_percent_id',
               'scaled_local_symmetric_percent_id',
               'query_align_len',
               'query_align_cov',
               'subject_align_len',
               'subject_align_cov',
               'bit_score']

    data_dict = df.to_dict('records')

    # Initialize PairWiseAligner with required parameters for model
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = mode

    final_data = []

    # Iterate through all pairs and calculate best alignment and metrics
    for row in data_dict:
        query = row['query']
        subject = row['subject']

        alignment = aligner.align(subject, query)
        best_alignment = max(alignment, key=lambda x: x.score)

        seq1_aligned = best_alignment[0]
        seq2_aligned = best_alignment[1]

        # Coverage
        subject_cov = sum(c != '-' for c in seq1_aligned) / len(seq1_aligned)
        query_cov = sum(c != '-' for c in seq2_aligned) / len(seq2_aligned)

        # Sequence length
        l_subject = len(query)
        l_query = len(subject)

        # Percent ids
        n_matches, n_gaps, n_columns, n_comp_gaps = get_matches_gaps(seq2_aligned, seq1_aligned)
        gap_comp_pct_id = gap_compressed_percent_id(n_matches, n_gaps, n_columns, n_comp_gaps)
        scaled_local_symmetric_percent_id = 2*n_matches / (l_subject + l_query)
        scaled_local_query_percent_id = n_matches / l_query

        # Collect calculated metrics into a list for addition to final_data
        new_row = [gap_comp_pct_id, scaled_local_query_percent_id,
                   scaled_local_symmetric_percent_id, l_query,
                   query_cov, l_subject, subject_cov, best_alignment.score,
                   row['qid'], row['sid']]
        final_data.append(new_row)

    # Construct temporary dataframe from collected metrics
    columns = metrics + ['query_id', 'subject_id']
    final_df = pd.DataFrame(final_data, columns=columns)

    # Merge final_df with sequences and ids from input df
    blast_df = df.merge(final_df, left_on=['qid', 'sid'],
                        right_on=['query_id', 'subject_id'])

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
