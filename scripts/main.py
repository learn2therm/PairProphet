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

from FAFSA.evaluate_model import evaluate_model
from FAFSA.train_val_wrapper import train_val_wrapper
from FAFSA.train_val_input_cleaning import columns_to_keep
from FAFSA.preprocessing import connect_db, build_fafsa
from FAFSA.structure import download_structures, run_fatcat


def model_dev(dbpath):

    #Ryans Component
    """
    Input: learn2thermDB (queuing format)
    Output: SQL db after filtering, and then each component can retrieve the data as needed.
    Params: Base environment + DuckDB
    """
    con, _ = connect_db(dbpath)
    build_fafsa(con)

    #Humood/Amin Component
        ##Humood:
    """
    Input: SQL db from Ryan, then retrieve only 1000 pairs
    Output: SQL db with appended Boolean from running HMMR remotely
    Params: Base environment + PyHMMER + Joblib
    """
    df = con.execute("""SELECT m_protein_seq, t_protein_seq, prot_pair_index, meso_pid, thermo_pid, meso_pdb, thermo_pdb""").df()
    test = con.execute("""SELECT pid, protein_seq FROM fafsa_proteins""")
    run_hmmr(df)
    parse_hmmer(path_to_csv)

    #Chau Component
    """
    Input: SQL db from Ryan, retrieve only PDB IDs and Uniprot IDs (sequences are not required for this component)
    Output: SQL db appended with Boolean on whether the structures are similar
    Params: Base environment + FATCAT
    """
    download_structures(df)
    run_fatcat(df)

    #Logan Component
    """
    Input: SQL db from Humood or Amin (with HMMR Boolean + prot_pair_ind + Jaccard score appended)
    Output: csv files with 6 columns (seq1, seq2, protein_match (Humood + Amin), protein_match(ML), Jaccard Score)
    Params: Base environment + iFeatureOmega dependencies 
    """
    #need to save model to use in next script
    model_output = train_val_wrapper(dataframe, feature_list=[])
    con.commit()
    con.close()
    return model_output

trained_model = model_dev(dbpath)[2]


def main(test_sequences):

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

    evaluate_model(output_path:str, trained_model, dataframe, target=[])

    return evaluation