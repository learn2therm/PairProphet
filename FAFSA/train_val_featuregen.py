"""
This module utilizes iFeatureOmega, a feature generation
package for proteins and nucleic acids.
"""

import iFeatureOmegaCLI
import pandas as pd
import Bio.SeqIO
import io
from io import StringIO

"""
Define a feature list from iFeatureOmega's catalog of descriptors.
Features have individual runtime <40 seconds when creating new dataframe.
Omitted some of the very high dimensional features (>1000-D).
"""

feature_list = ['AAC', 'GAAC', 'DistancePair',
                'DPC type 1', 'CTDC', 'CTDT', 'CTDD', 'CTriad',
                'CKSAAGP type 1', 'PseKRAAC type 1', 'APAAC', 'QSOrder']


def get_fasta_from_dataframe(dataframe, output_file: str):
    # adjust this to write function with BioPython
    """
    Function converts sequences and ID's from
    pandas dataframe into .fasta file that can
    be read by iFeatureOmega package.

    Parameters
    ----------
    dataframe : Pandas dataframe

    Returns
    ----------
    Writes .fasta file to current directory.
    """

    with open(output_file, 'w') as f:
        for _, row in dataframe.iterrows():
            f.write(
                '>{}\n{}\n'.format(
                    (row['meso_index']),
                    row['m_protein_seq']))
    return output_file


def get_protein_descriptors(fasta: str, descriptors=[]):
    """
    Generates features from a protein sequence

    Parameters
    ----------
    Fasta file with protein sequences.

    Returns
    -------
    Vector of descriptors
    """

    # create iProtein object
    protein = iFeatureOmegaCLI.iProtein(fasta)

    # not sure why we need this yet. Right now it is stored in local directory.
    params = protein.import_parameters('protein_parameters.json')

    protein_descriptors = {}

    for descriptor in descriptors:
        protein.get_descriptor(descriptor)
        protein_descriptors.update({f'{descriptor}': protein.encodings})

    return protein_descriptors


def create_new_dataframe(dataframe, output_file, descriptors=[]):
    """
    Creates new dataframe with descriptors added.

    Parameters
    ----------
    Pandas dataframe, list of descriptors as strings, output file name.

    Returns
    -------
    Dataframe including vector(s) of descriptors
    """

    fasta = get_fasta_from_dataframe(dataframe, output_file)

    feature_dict = get_protein_descriptors(fasta, descriptors)

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

    return df


"""
Functions below only necessary if input is fasta file.
"""


def remove_fasta_description(filename: str):
    """
    Removes description from fasta file so that iProtein can read the input.
    """

    # assign unwanted string to object
    string_to_remove = "<unknown description>"

    # open file
    with open(filename, "r") as file:
        content = file.read()

    # Remove the string
    new_content = content.replace(string_to_remove, "")

    # overwrite file without string
    with open(filename, "w") as file:
        seq = file.write(new_content)

    return seq


"""
Function below is necessary for inputs with
two protein sequences.
"""

def get_fasta_from_dataframe(dataframe, output_file_a, output_file_b):
    """
    adjust this to write function with BioPython
    separate functions for each of the input sequences
    in training, seq_a = meso and seq_b = thermo
    """
    
    
    #meso sequence to fasta
    with open(output_file_a, 'w') as f:
        for _, row in df.iterrows():
            f.write('>{}\n{}\n'.format((row['prot_pair_index']), row['m_protein_seq']))
    
    #thermo sequence to fasta
    with open(output_file_b, 'w') as f:
        for _, row in df.iterrows():
            f.write('>{}\n{}\n'.format((row['prot_pair_index']), (row['t_protein_seq'])))
   
    #return output files
    return [output_file_a, output_file_b]

def create_new_dataframe(dataframe, output_files, descriptors=[]):
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

    def compute_descriptor_difference(fasta_files, descriptors=[]):
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
            feature_dict[key] = desc_a[key] - desc_b[key]

        return feature_dict

    feature_dict = compute_descriptor_difference(fasta_files, descriptors)

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

    return df
