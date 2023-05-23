"""
This module utilizes iFeatureOmega, a feature generation
package for proteins and nucleic acids.
"""

import iFeatureOmegaCLI
import pandas as pd
# import Bio.SeqIO
# import io
# from io import StringIO


feature_list = ['AAC', 'GAAC', 'DistancePair',
                'CTDC', 'CTDT', 'CTDD', 'CTriad', 'GDPC type 1', 'GDPC type 2',
                'CKSAAGP type 1', 'CKSAAGP type 2', 'PseKRAAC type 2', 'PseKRAAC type 3A',
                'PseKRAAC type 7', 'PseKRAAC type 9', 'Geary', 'APAAC', 'QSOrder']


def get_fasta_from_dataframe(
        dataframe, output_file_a: str, output_file_b: str):
    """
    Generates fasta file type from pandas dataframe.

    Parameters
    ----------
    -Dataframe (pandas dataframe)
    -Names of output fasta files (str)

    Returns
    -------
    Two fasta files with protein sequences and prot_pair_index

    Note: Want to adjust this to write function with BioPython
    separate functions for each of the input sequences
    in training, seq_a = meso and seq_b = thermo
    """
    # meso sequence to fasta
    with open(output_file_a, 'w') as f:
        for _, row in dataframe.iterrows():
            f.write(
                '>{}\n{}\n'.format(
                    (row['prot_pair_index']),
                    row['m_protein_seq']))

    # thermo sequence to fasta
    with open(output_file_b, 'w') as f:
        for _, row in dataframe.iterrows():
            f.write(
                '>{}\n{}\n'.format(
                    (row['prot_pair_index']),
                    (row['t_protein_seq'])))

    # return output files
    return [output_file_a, output_file_b]


def get_protein_descriptors(fasta_file: str, descriptors=[]):
    """
    Generates features from a protein sequence

    Parameters
    ----------
    Fasta file with amino acid sequences.

    Returns
    -------
    Vector of descriptors
    """

    # create iProtein object
    protein = iFeatureOmegaCLI.iProtein(fasta_file)

    # not sure why we need this yet. Right now it is stored in local directory.
    params = protein.import_parameters('protein_parameters.json')

    protein_descriptors = {}

    for descriptor in descriptors:
        protein.get_descriptor(descriptor)
        protein_descriptors.update({f'{descriptor}': protein.encodings})

    # make sure output is a dictionary of correct length
    assert "dict" in str(type(protein_descriptors))
    assert len(protein_descriptors) == len(descriptors)

    return protein_descriptors


def clean_new_dataframe(dataframe):
    """Creates two new columns of bit score
    normalized by the protein length.

    Returns
    -------
    Pandas dataframe
    """

    # drop indexing columns created by feature gen
    dataframe = dataframe.drop(
        columns=dataframe.columns[dataframe.columns.str.contains('index|Unnamed')])

    assert (dataframe.filter(like='index|Unnamed').shape)[1] == 0

    # turn inf into NaN
    dataframe = dataframe.replace([np.inf, -np.inf], np.nan)

    # assert NaN's are removed
    nan_counts = df.isna().sum()
    assert nan_counts.unique() == [0]

    return dataframe


def create_new_dataframe(dataframe, output_files: list, descriptors=[]):
    """
    Creates new dataframe with descriptors added.

    Parameters
    ----------
    Pandas dataframe, list of descriptors as strings, output file name.

    Returns
    -------
    Dataframe including vector(s) of descriptors
    """

    fasta_files = get_fasta_from_dataframe(
        dataframe, output_files[0], output_files[1])

    def compute_descriptor_ratio(fasta_files, descriptors=[]):
        """
        Generates dictionary of descriptors for each of the two input sequences.
        Computes the difference between each instance of a descriptor.

        Parameters
        ----------
        List of two fasta files (str) and list of descriptors (str).

        Returns
        -------
        Dictionary with difference between descriptors for each of the
        input sequences.
        """
        desc_a = get_protein_descriptors(fasta_files[0], descriptors)
        desc_b = get_protein_descriptors(fasta_files[1], descriptors)

        feature_dict = {}

        for key in desc_a:

            if 'AAC' in key:
                feature_dict[key] = desc_a[key] - desc_b[key]
            elif 'GAAC' in key:
                feature_dict[key] = desc_a[key] - desc_b[key]
            else:
                feature_dict[key] = desc_a[key] / desc_b[key]

        return feature_dict

    feature_dict = compute_descriptor_ratio(fasta_files, descriptors)

    df = dataframe.reset_index()

    for desc in descriptors:

        feature_dict[desc].index = feature_dict[desc].index.astype(int)
        features = feature_dict[desc].reset_index()

        df = pd.merge(
            df,
            features,
            how='outer',
            left_index=True,
            right_index=True)

    df = clean_new_dataframe(df)

    return df
