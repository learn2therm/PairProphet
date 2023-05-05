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

    #Humood/Amin Component

    #Chau Component
    """
    Input: Dataframe with PDB IDs and Uniprot IDs (sequences are not required for this component)
    Output: Boolean on whether the structures are similar
    Params: Base environment + BioPython + FATCAT
    """
    #Logan Component
    """
    Input: Dataframe with BLAST metrics/sequences from Ryan (pandas dataframe) 
    + Boolean generated from Humood (.csv)
    Output: Pandas dataframe with 4 columns (seq1, seq2, protein_match, Jaccard Score)
    Params: Base environment + iFeatureOmega dependencies 
    """
    model_output = rf_wrapper(blast_seq, boolean1, boolean2)

    return model_output


def main():

    """
    Function for user to interact with.
    """