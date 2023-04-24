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

# def remove_fasta_description(filename:str):

#     """
#     Removes description from fasta file so that iProtein can read the input.
#     """

#     #assign unwanted string to object
#     string_to_remove = "<unknown description>"

#     #open file
#     with open(filename, "r") as file:
#         content = file.read()

#     # Remove the string
#     new_content = content.replace(string_to_remove, "")

#     #overwrite file without string
#     with open(filename, "w") as file:
#         seq = file.write(new_content)

#     return seq

# def fasta_to_descriptors(fasta:str, descriptors=[]):

#     """This function is for generating descriptors from a fasta file.
#     Not necessary if sequences come from dataframe.

#     Parameters
#     ----------
#     Fasta file with protein sequences.

#     Returns
#     -------


#     #remove description from fasta file
#     remove_fasta_description(fasta)

#     #return protein descriptors
#     return get_protein_descriptors(fasta, descriptors=descriptors)
