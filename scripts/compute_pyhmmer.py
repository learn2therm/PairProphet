'''
The script executes HMMER against pfam and generate output files

The packages you need to run this script are the following:
- biopython
- pyHMMER (https://pyhmmer.readthedocs.io/en/stable/)
- pandas
- joblib

You also need to have:
- pfam db locally
- protein db
'''
# system dependecies
import os
from pathlib import Path
import sys
from typing import Union

# job stuff
import logging
from joblib import delayed, Parallel


# library dependencies
import pandas as pd
from tqdm import tqdm

# biopython
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# pyhmmer
import pyhmmer


# Paths (used to need as I was planning to use this script for hyak)
PFAM_PATH = Path("/Users/humoodalanzi/pfam/Pfam-A.hmm")  # ./Pfam-A.hmm
# ./learn2therm_sample_50k.csv
SAMPLE_DB_PATH = Path("../notebooks/learn2therm_sample_50k_exploration.csv")



def hmmpress_hmms(hmms_path, pfam_data_folder):
    """
    Presses the HMMs in the given HMM database and stores the resulting files in a specified directory.

    Parameters
    ----------
    hmmdb_path : str
        Path to the HMM database.
    pfam_data_folder : str, optional
        Path to the directory where the HMMs should be stored.

    Returns
    -------
    None

    Notes
    -----
    This function uses HMMER's hmmpress program to compress the HMMs in the given HMM database and
    stores the resulting files in the specified directory for faster access during future HMMER runs.
    If the specified directory does not exist, it will be created.
    """
    hmms = pyhmmer.plan7.HMMFile(hmms_path)
    pyhmmer.hmmer.hmmpress(hmms, pfam_data_folder)


def prefetch_targets(hmms_path: str):
    """
    Prefetch HMM profiles from a given HMM database.

    Parameters
    ----------
    hmms_path : str
        Path to the HMM database.

    Returns
    -------
    targets : pyhmmer.plan7.OptimizedProfileBlock
        The HMM profiles loaded from the database.
    """
    # amino acid alphabet and prefetched inputs
    amino_acids = pyhmmer.easel.Alphabet.amino()
    optimized_profiles = list(pyhmmer.plan7.HMMPressedFile(hmms_path))
    targets = pyhmmer.plan7.OptimizedProfileBlock(
        amino_acids, optimized_profiles)  
    return targets


def save_sequences_to_fasta(sequences: pd.core.frame.DataFrame, inputname: str = "input"):
    """
    Returns a list of SeqRecord objects and creates a corresponding input Fasta of them

    Parameters:
    ------------
    sequences : pandas.core.frame.DataFrame
        a dataframe with string amino acid sequences in a 'protein_seq' column
    input name : str, default = 'input'
        a name for the input fasta file


    Returns:
    ------------
    file : TextIOWrapper
        the input fasta file created from the list of SeqRecord objects

    Raises
    -------
    ValueError :
        if the input dataframe is empty
    AttributeError :
        if any of the sequences are invalid
    """
    # ensure input file has .fasta extension
    if not inputname.endswith('.fasta'):
        inputname = f"{os.path.splitext(inputname)[0]}.fasta"

    # check if input is empty
    if sequences.empty:
        raise ValueError("Input dataframe is empty")

    # check if sequences are valid
    for seq in sequences['protein_seq']:
        try:
            Seq(seq)
        except BaseException as exc:
            raise AttributeError("Invalid sequence") from exc

    # function
    records = []
    for index, seq in sequences.itertuples():
        try:
            record = SeqRecord(Seq(seq), id=str(index))
            records.append(record)
        except AttributeError as exc:
            raise AttributeError(f"Invalid sequence: {seq}") from exc

    # raise error if seq not valid
    if not records:
        raise AttributeError("No valid sequences found in input")

    with open(inputname, "w", encoding="utf-8") as file:
        SeqIO.write(records, file, "fasta")
    return file

def run_pyhmmer(
        input_file: str,
        hmms_path: str,
        prefetch: bool=False,
        output_file: str = None,
        cpu: int = 4,
        eval_con: float = 1e-10):
    """
    Run HMMER's hmmscan program on a set of input sequences using with HMMs from a database.

    Parameters
    ----------
    input_file : str
        Path to the input sequence file.
    hmms_path : str
        Path to the HMM database.
    prefetch : bool, optional
        Specifies how the HMM are stored in meomry.
    output_file : str, optional
        Path to the output file if the users wants to write the file.
    cpu : int, optional
        The number of CPUs to use. Default is 4.
    eval_con : float, optional
        E-value threshold for domain reporting. Default is 1e-10.

    Returns
    -------
    all_hits : pyhmmer.plan7.TopHits or domtblout file
        If the output_file has a name, it will be written to a domtblout file.
        Otherwise, the user will get a list of pyhmmeer TopHits objects.

    Notes
    -----
    This function runs HMMER's hmmscan program on a set of input sequences
    using HMMs from a given database.
    The function supports two modes: normal mode and prefetching mode.
    In normal mode, the HMMs are pressed and stored in a directory before execution.
    In prefetching mode, the HMMs are kept in memory for faster search.
    """ 
    # ensure input file has .fasta extension
    if not input_file.endswith('.fasta'):
        input_file = f"{os.path.splitext(input_file)[0]}.fasta"

    # ensure output_file has .domtblout extension
    if not output_file.endswith('.domtblout'):
        output_file = f"{os.path.splitext(output_file)[0]}.domtblout"

    # HMM profile modes
    if prefetch:
        targets = prefetch_targets(hmms_path)
    else:
        targets = pyhmmer.plan7.HMMFile("../data/pfam/.h3m")
    
    # HMMscan execution with or without save_out
    with pyhmmer.easel.SequenceFile(input_file, digital=True) as seqs:
        all_hits = list(pyhmmer.hmmer.hmmscan(seqs, targets, cpus=cpu, E=eval_con))
        # check if we should save the output
        if output_file is not None:
            with open(output_file, "wb") as dst:
                for i, hits in enumerate(
                    pyhmmer.hmmer.hmmscan(
                        seqs, targets, cpus=cpu, E=eval_con)):
                    hits.write(dst, format="domains", header=i == 0)
    return all_hits


def parse_pyhmmer(all_hits):
    """
    Parses the TopHit pyhmmer object getting the query and accession IDs and saves to a DataFrame

    Parameters
    ----------
    all_hits : list
        A list of TopHit objects from pyhmmer.

    Returns
    -------
    pandas.DataFrame
        A dataframe containing the query and accession IDs.
    """
    # initialize an empty list to store the data
    data = []
    
    # iterate over each protein hit
    for top_hits in all_hits:
        for hit in top_hits:
            # extract the query and accession IDs and decode the query ID
            query_id = hit.hits.query_name.decode('utf-8')
            accession_id = hit.accession.decode('utf-8')
            
            # append the data to the list
            data.append([query_id, accession_id])
    
    # create the DataFrame from the list
    df = pd.DataFrame(data, columns=["query_id", "accession_id"])
    
    # group the accession IDs by query ID and join them into a single string separated by ";"
    df = df.groupby("query_id")["accession_id"].apply(lambda x: ";".join(x)).reset_index()
    
    return df


if __name__ == '__main__':
    # Logging
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    fh = logging.FileHandler('hmmer.log', mode='w')
    formatter = logging.Formatter(
        '%(filename)-12s %(asctime)s;%(funcName)-12s: %(levelname)-8s %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    logger.info("TEST LOG")

    # Set up parallel processing and parsing
    total_size = int(sys.argv[1])  # Number of total sequences read
    chunk_size = int(sys.argv[2])  # Number of sequences to process in each chunk
    njobs = int(sys.argv[3])  # Number of parallel processes to use

    logger.info('Parallel processing parameters obtained')

    # Data prep and processing

    # reading the data
    # sample
    df_sample = pd.read_csv(SAMPLE_DB_PATH, index_col=0)
    logger.info('Loaded database')

    # separating the meso and thermo
    meso_seq_db = df_sample[["meso_index", "m_protein_seq"]]
    thermo_seq_db = df_sample[["thermo_index", "t_protein_seq"]]
    logger.info('Data seperated into t and m')

    # processing meso to be suitable for HMMER
    meso_seq_list = meso_seq_db.set_index("meso_index").iloc[:total_size]
    meso_seq_list.index.name = None
    meso_seq_list.rename({'m_protein_seq': 'protein_seq'},
                         axis="columns", inplace=True)

    # processing thermo to be suitable for HMMER
    thermo_seq_list = thermo_seq_db.set_index("thermo_index").iloc[:total_size]
    thermo_seq_list.index.name = None
    thermo_seq_list.rename(
        {'t_protein_seq': 'protein_seq'}, axis="columns", inplace=True)

    logger.info('Sampled t and m data')

    # press the HMM db
    hmmpress_hmms(PFAM_PATH, "../data/pfam/")

    logger.info('Pressed HMM DB')

    # create worker function to scatter 
    def worker_function(chunk_index, sequences, which):
        """
        A wrapping function that runs and parses pyhmmer in chunks

        Parameters
        ----------
        chunk_index : int
            number of sequences chunks
        sequences : str
            a list of dataframe containing protein sequences
        which : bool
            class of proteins
        """
        # define paths for input and output files
        input_file_path = f"./results/{which}_input_{chunk_index}"
        output_file_path = f"./results/{which}_output_{chunk_index}"

        # convert sequences to FASTA files
        save_sequences_to_fasta(sequences, input_file_path)

        hits = run_pyhmmer(
            input_file=input_file_path,
            hmms_path=PFAM_PATH,
            prefetch=True,
            output_file=output_file_path,
            cpu=1,
            eval_con=1e-5)

        # Parse pyhmmer output and save to CSV file
        accessions_parsed = parse_pyhmmer(all_hits=hits)
        accessions_parsed.to_csv(f"./results/{which}_result_{chunk_index}.csv", index=False)

    # chunking the data to chunk_size sequence bits (change if sample or all
    # proteins)
    meso_chunks = [meso_seq_list[i:i + chunk_size]
                   for i in range(0, len(meso_seq_list), chunk_size)]
    thermo_chunks = [thermo_seq_list[i:i + chunk_size]
                     for i in range(0, len(thermo_seq_list), chunk_size)]

    logger.info('Sequence chunking done')

    # parallel computing on how many CPUs (n_jobs=)
    logger.info('Running pyhmmer in parallel on all chunks')

    Parallel(n_jobs=njobs)(delayed(worker_function)(chunk_index, sequences , "meso")
                               for chunk_index, sequences  in enumerate(meso_chunks))
    Parallel(n_jobs=njobs)(delayed(worker_function)(chunk_index, sequences , "thermo")
                               for chunk_index, sequences in enumerate(thermo_chunks))

    logger.info('Parallelization complete')
