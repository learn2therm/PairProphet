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

def hmmpress_hmms(hmms_path, hmm_db_path):
    """
    Presses the HMMs in the given HMM database and stores the resulting files in a specified directory.
    Parameters
    ----------
    hmmdb_path : str
        Path to the HMM database.
    hmm_db_path : str, optional
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
    pyhmmer.hmmer.hmmpress(hmms, hmm_db_path)

def prefetch_targets(hmm_db_path: str):
    """
    Prefetch HMM profiles from a given HMM database.
    
    Parameters
    ----------
    hmmdb : str
        Path to the HMM database.
    Returns
    -------
    targets : pyhmmer.plan7.OptimizedProfileBlock
        The HMM profiles loaded from the database.

    """
    # amino acid alphabet and prefetched inputs
    aa = pyhmmer.easel.Alphabet.amino()
    optimized_profiles = list(pyhmmer.plan7.HMMPressedFile(hmm_db_path))
    targets = pyhmmer.plan7.OptimizedProfileBlock(
        aa, optimized_profiles) 
    return targets

def run_pyhmmer(
        input_file: str,
        hmm_db_path: str,
        prefetch: bool = False,
        output_file: str = None,
        cpu: int = 4,
        eval_con: float = 1e-10):
    """
    Runs HMMER's hmmscan program on a set of input sequences using HMMs from a given database.

    Parameters
    ----------
    input_file : str
        Path to the input sequence file.
    hmm_db_path : str
        Path to the HMM database.
    prefetch: bool,
    output_file : str
        Path to the output file.
    cpu : int, optional
        The number of CPUs to use. Default is 4.
    save_out : bool, optional
        Whether to save the output to file. Default is False.
    eval_con : float, optional
        E-value threshold for domain reporting. Default is 1e-10.

    Returns
    -------
    pyhmmer.plan7.TopHits or None
        The output hits if `save_out` is False, otherwise None.

    Raises
    ------
    ValueError
        If the input dataframe is empty.
    AttributeError
        If any of the sequences are invalid.

    Notes
    -----
    This function runs HMMER's hmmscan program on a set of input sequences
    using HMMs from a given database.
    The function supports two modes: normal mode and prefetching mode.
    In normal mode, the HMMs are pressed and stored in a directory before execution.
    In prefetching mode, the HMMs are kept in memory for faster search.
    """  
    # Ensure output_file has .domtblout extension
    if not output_file.endswith('.domtblout'):
        output_file = f"{os.path.splitext(output_file)[0]}.domtblout"

    if prefetch:
        targets = prefetch_targets(hmm_db_path)
    else:
        targets = pyhmmer.plan7.HMMPressedFile(hmm_db_path)

    # HMMscan execution with or without save_out
    with pyhmmer.easel.SequenceFile(input_file, digital=True) as seqs:
        all_hits = list(pyhmmer.hmmer.hmmscan(seqs, targets, cpus=cpu, E=eval_con))
        # check if we should save to output
        if output_file is not None:
            with open(output_file, "wb") as dst:
                for i, hits in enumerate(all_hits):
                    hits.write(dst, format="domains", header= i==0)
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

def save_sequences_to_fasta(sequences: pd.core.frame.DataFrame, inputname: str = "input"):
    """
    Returns a list of SeqRecord objects and creates a corresponding input Fasta of them
    Parameters:
    ------------
    list : pandas.core.frame.DataFrame
        a dataframe with string amino acid sequences in a 'seq' column
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
    # make sure it ends with fasta
    if not input_file.endswith('.fasta'):
        input_file = f"{os.path.splitext(input_file)[0]}.fasta"

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

    # press the db
    hmmpress_hmms(PFAM_PATH, "../data/pfam/")

    # create worker function to scatter
    def worker_function(sequences, chunk_index, which):
        input_file_path = f"./results/{which}_input_{chunk_index}.fasta"
        output_file_path = f"./results/{which}_output_{chunk_index}.domtblout"

        save_sequences_to_fasta(sequences, input_file_path)

        hits = run_pyhmmer(
            input_file_path,
            "../data/pfam/pfam.h3m", # not sure if pathing is correct here, need to check paths since this was designed for your computer
            output_file=output_file_path,
            prefetch=True,
            cpu=1,
            eval_con=1e-5,
        )

        accessions_parsed = parse_pyhmmer(hits)
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

    Parallel(n_jobs=njobs)(delayed(worker_function)(sequences, chunk_index, "meso")
                               for sequences, chunk_index in enumerate(meso_chunks))
    Parallel(n_jobs=njobs)(delayed(worker_function)(sequences, chunk_index, "thermo")
                               for sequences, chunk_index in enumerate(thermo_chunks))

    logger.info('Parallelization complete')
