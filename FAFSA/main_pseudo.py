"""
Main script wrapper.
This is pseudocode for our parent script.
Directory: ValidProt.
Input: 
Output: 
Note: We have two wrapping functions, but 
these can be developed as individual scripts.
Runtime: 
"""

def model_dev(sequences):

    #Ryans Component
    """
    Input: learn2thermDB (queuing format)
    Output: SQL db after filtering, and then each component can retrieve the data as needed.
    Params: Base environment + DuckDB
    """

    #Humood/Amin Component
        ##Humood:
    """
    Input:
    Output:
    Params:
    """
        ##Amin:
    """
    Input: SQL db from Ryan, then retrieve only 1000 pairs
    Output: SQL db with appended Boolean from running HMMR remotely
    Params: Base environment + HMMR API
    """
    #Chau Component
    """
    Input: SQL db from Ryan, retrieve only PDB IDs and Uniprot IDs (sequences are not required for this component)
    Output: SQL db appended with Boolean on whether the structures are similar
    Params: Base environment + BioPython + FATCAT
    """
    #Logan Component
    """
    Input: SQL db from Humood or Amin (with HMMR Boolean + prot_pair_ind + Jaccard score appended)
    Output: csv files with 6 columns (seq1, seq2, protein_match (Humood + Amin), protein_match(ML), Jaccard Score)
    Params: Base environment + iFeatureOmega dependencies 
    """
    model_output = rf_wrapper(blast_seq, boolean1, boolean2)

    return model_output


def main():

    """
    Function for user to interact with.
    """
    #Ryans Component
    """
    Input: csv file (list of sequences) + pairing and PDBids and UniprotIDs are optional
    Output:chunked csv for HMMR and structure components.
        Note: csv with specific parameters will be generated for specific component.
    Params:
    """
    #Logan Component
    """
    Input: SQL db from Humood or Amin + Chau (with HMMR Boolean + prot_pair_ind + Jaccard score appended + FATCAT)
    Output: csv files with 6 columns (seq1, seq2, protein_match (Humood + Amin), protein_match (Chau) protein_match(ML), Jaccard Score)
    Params: Base environment + iFeatureOmega dependencies 
    """