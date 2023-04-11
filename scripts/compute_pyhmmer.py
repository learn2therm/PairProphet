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


# job stuff
import logging
from joblib import delayed, Parallel


# library dependencies
import pandas as pd


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


def run_hmmer(
        seqs: pd.core.frame.DataFrame,
        input_file: str,
        hmm: str,
        output_file: str,
        cpu: int = 4,
        prefetching=False,
        save_out=False,
        eval_con: float = 1e-10):
    """
    Runs HMMER's hmmscan program on a set of input sequences using HMMs from a given database.

    Parameters
    ----------
    seqs : pandas.core.frame.DataFrame
        A dataframe with string amino acid sequences in a 'seq' column.
    input_file : str
        Path to the input sequence file.
    hmm : str
        Path to the HMM database.
    output_file : str
        Path to the output file.
    cpu : int, optional
        The number of CPUs to use. Default is 4.
    prefetching : bool, optional
        Whether to use prefetching for faster search. Default is False.
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
    # generate meso and thermo files
    read_seq(seqs, input_file)

    # place files into HMMER/pfam
    run_pyhmmer(
        hmm,
        input_file,
        output_file,
        cpu,
        prefetching,
        save_out,
        eval_con)


def read_seq(lists: pd.core.frame.DataFrame, inputname: str = "input"):
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
    # check if input is empty
    if lists.empty:
        raise ValueError("Input dataframe is empty")

    # check if sequences are valid
    for seq in lists['protein_seq']:
        try:
            Seq(seq)
        except BaseException as exc:
            raise AttributeError("Invalid sequence") from exc

    # function
    records = []
    for index, seq in lists.itertuples():
        try:
            record = SeqRecord(Seq(seq), id=str(index))
            records.append(record)
        except AttributeError as exc:
            raise AttributeError(f"Invalid sequence: {seq}") from exc

    # raise error if seq not valid
    if not records:
        raise AttributeError("No valid sequences found in input")

    with open(f"{inputname}.fasta", "w", encoding="utf-8") as file:
        SeqIO.write(records, file, "fasta")
    return file


def run_pyhmmer(
        hmmdb: str,
        input_file: str,
        output_file: str,
        cpu: int = 4,
        prefetching=False,
        save_out=False,
        eval_con: float = 1e-10):
    """
    Run hmmscan on input sequences with HMMs from a database.

    Parameters
    ----------
    hmmdb : str
        Path to the HMM database.
    input_file : str
        Path to the input sequence file.
    output_file : str
        Path to the output file.
    cpu : int, optional
        The number of CPUs to use. Default is 4.
    prefetching : bool, optional
        Whether to use prefetching for faster search. Default is False.
    save_out : bool, optional
        Whether to save the output to file. Default is False.
    eval_con : float, optional
        E-value threshold for domain reporting. Default is 1e-10.

    Returns
    -------
    all_hits : pyhmmer.plan7.TopHits or None
        The output hits if `save_out` is False, otherwise None.

    Notes
    -----
    This function runs HMMER's hmmscan program on a set of input sequences
    using HMMs from a given database.
    The function supports two modes: normal mode and prefetching mode.
    In normal mode, the HMMs are pressed and stored in a directory before execution.
    In prefetching mode, the HMMs are kept in memory for faster search.
    """
    # Create hmms
    hmms = pyhmmer.plan7.HMMFile(hmmdb)
    # press hmms and store them in the pfam data folder or w/e destination
    if not os.path.exists(
        os.path.join(
            "../data/pfam/",
            os.path.basename(hmmdb) +
            ".h3m")):
        pyhmmer.hmmer.hmmpress(hmms, "../data/pfam/")

    # Ensure input_file has .fasta extension
    if not input_file.endswith('.fasta'):
        input_file = f"{os.path.splitext(input_file)[0]}.fasta"
    # Ensure output_file has .domtblout extension
    if not output_file.endswith('.domtblout'):
        output_file = f"{os.path.splitext(output_file)[0]}.domtblout"

    # amino acid alphabet and prefetched inputs
    aa = pyhmmer.easel.Alphabet.amino()
    optimized_profiles = list(pyhmmer.plan7.HMMPressedFile(hmmdb))
    targets = pyhmmer.plan7.OptimizedProfileBlock(
        aa, optimized_profiles) if prefetching else pyhmmer.plan7.HMMFile("../data/pfam/.h3m")

    # HMMscan execution with or without save_out
    with pyhmmer.easel.SequenceFile(input_file, digital=True) as seqs:
        if save_out:
            with open(output_file, "wb") as dst:
                for i, hits in enumerate(
                    pyhmmer.hmmer.hmmscan(
                        seqs, targets, cpus=cpu, E=eval_con)):
                    hits.write(dst, format="domains", header=i == 0)
        else:
            all_hits = pyhmmer.hmmscan(seqs, targets, cpus=cpu, E=eval_con)

    return all_hits if not save_out else None


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
    meso_seq_list = meso_seq_db.set_index("meso_index").iloc[:500]
    meso_seq_list.index.name = None
    meso_seq_list.rename({'m_protein_seq': 'protein_seq'},
                         axis="columns", inplace=True)

    # processing thermo to be suitable for HMMER
    thermo_seq_list = thermo_seq_db.set_index("thermo_index").iloc[:500]
    thermo_seq_list.index.name = None
    thermo_seq_list.rename(
        {'t_protein_seq': 'protein_seq'}, axis="columns", inplace=True)

    logger.info('Sampled t and m data')

    # Execution
    def run_hmmer_parallel(chunk_index, seq, which):
        """Executes run_hmmer in parrallel"""
        run_hmmer(
            seqs=seq,
            input_file=f"{which}_input_{chunk_index}",
            hmm=PFAM_PATH,
            output_file=f"{which}_output_{chunk_index}.domtblout",
            cpu=1,
            prefetching=True,
            save_out=True
        )

    # chunking the data to 100 sequence bits (change if sample or all proteins)
    meso_chunks = [meso_seq_list[i:i + 100]
                   for i in range(0, len(meso_seq_list), 100)]
    thermo_chunks = [thermo_seq_list[i:i + 100]
                     for i in range(0, len(thermo_seq_list), 100)]

    # parallel computing on how many CPUs (n_jobs=)
    Parallel(n_jobs=1)(delayed(run_hmmer_parallel)(i, chunk, "meso")
                       for i, chunk in enumerate(meso_chunks))
    Parallel(n_jobs=1)(delayed(run_hmmer_parallel)(i, chunk, "thermo")
                       for i, chunk in enumerate(thermo_chunks))
