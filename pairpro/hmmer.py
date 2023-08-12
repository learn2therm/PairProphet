"""
The following are various importable code for running HMMER
either locally via pyhmmer or via Interpro's API

The local version runs in parrallel using joblib.
The API version runs in parrallel using asyncio.

The local version is faster, but the API version is more accessible.

TODO:
    - Add error handling/tests
    - Update parsing function documentation
"""
# system dependecies
import asyncio
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import csv
import glob
import json
import logging
import math
import nest_asyncio
import os
import requests
import time
import urllib.parse
import tempfile
from typing import Dict, List, Tuple, Union

# library dependencies
import duckdb as ddb
import httpx
from joblib import Parallel, delayed
import pandas as pd
import pyhmmer

# local dependencies
import pairpro.utils

logger = logging.getLogger(__name__)


# Local HMMER

def hmmpress_hmms(hmms_path, pfam_data_folder):
    """
    Presses the HMMs in the given HMM database and stores the resulting files in a specified directory.

    Args:
        hmms_path (str): Path to the HMM database.
        pfam_data_folder (str): Path to the directory where the HMMs should be stored.

    Returns:
        None

    Notes:
        This function uses HMMER's hmmpress program to compress the HMMs in the given HMM database
        and stores the resulting files in the specified directory for faster access during future HMMER runs.
        If the specified directory does not exist, it will be created.
    """
    hmms = pyhmmer.plan7.HMMFile(hmms_path)
    pyhmmer.hmmer.hmmpress(hmms, pfam_data_folder)


def prefetch_targets(hmms_path: str):
    """
    Prefetch HMM profiles from a given HMM database.

    Args:
        hmms_path (str): Path to the pressed HMM database.

    Returns:
        targets (pyhmmer.plan7.OptimizedProfileBlock): The HMM profiles loaded from the database.
    """
    # amino acid alphabet and prefetched inputs
    amino_acids = pyhmmer.easel.Alphabet.amino()
    optimized_profiles = list(pyhmmer.plan7.HMMPressedFile(hmms_path))
    targets = pyhmmer.plan7.OptimizedProfileBlock(
        amino_acids, optimized_profiles)
    return targets


def save_to_digital_sequences(dataframe: pd.DataFrame):
    """
    Save protein sequences from a DataFrame to a digital sequence block.

    Args:
        dataframe (pd.DataFrame): DataFrame containing PIDs (Protein IDs) and sequences.

    Returns:
        pyhmmer.easel.DigitalSequenceBlock: A digital sequence block containing the converted sequences.
    """
    # Create empty list
    seqlist = []

    # Establish pyhmmer alphabet
    amino_acids = pyhmmer.easel.Alphabet.amino()

    # Convert proteins in dataframe to suitable format
    for _, row in dataframe.iterrows():
        pid = bytes(row['pid'], encoding='utf-8')
        seq_str = row['protein_seq']
        sequences = pyhmmer.easel.TextSequence(name=pid, sequence=seq_str)
        sequences = sequences.digitize(amino_acids)
        seqlist.append(sequences)

    # Convert so SequenceBlocks
    seqblock = pyhmmer.easel.DigitalSequenceBlock(amino_acids, seqlist)

    return seqblock


def run_pyhmmer(
        seqs: Union[pyhmmer.easel.DigitalSequenceBlock, str],
        hmms_path: str = None,
        pressed_path: str = None,
        prefetch: Union[bool, pyhmmer.plan7.OptimizedProfileBlock] = False,
        output_file: str = None,
        cpu: int = 4,
        scan: bool=True,
        eval_con: float = 1e-10,
        **kwargs
        ):
    """
    Run HMMER's hmmscan program on a set of input sequences using HMMs from a database.

    Args:
        seqs (pyhmmer.easel.DigitalSequenceBlock): Digital sequence block of input sequences.
        hmms_path (str): Path to the HMM database.
        pressed_path (str): Path to the pressed HMM database.
        prefetch (bool, optional): Specifies whether to use prefetching mode for HMM storage. Defaults to False.
        output_file (str, optional): Path to the output file if the user wants to write the file. Defaults to None.
        cpu (int, optional): The number of CPUs to use. Defaults to 4.
        scan (bool, optional): Specifies whether to run hmmscan or hmmsearch. Defaults to True.
        eval_con (float, optional): E-value threshold for domain reporting. Defaults to 1e-10.

    Returns:
        Union[pyhmmer.plan7.TopHits, str]: If output_file is specified, the function writes the results to a domtblout file and returns the file path.
        Otherwise, it returns a list of pyhmmer.plan7.TopHits objects.

    Notes:
        This function runs HMMER's hmmscan program on a set of input sequences using HMMs from a given database.
        The function supports two modes: normal mode and prefetching mode.
        In normal mode, the HMMs are pressed and stored in a directory before execution.
        In prefetching mode, the HMMs are kept in memory for faster search.
    """
    # check if both hmms_path and pressed_path are specified
    if not hmms_path and not pressed_path:
        raise ValueError("Must specifity one of hmm path (to a .hmm file) or pressed_path (containing .h3m, etc.)")
    
    # ensure output_file has .domtblout extension
    if output_file is not None and not output_file.endswith('.domtblout'):
        output_file = f"{os.path.splitext(output_file)[0]}.domtblout"

    # HMM profile modes
    if prefetch:
        if isinstance(prefetch, pyhmmer.plan7.OptimizedProfileBlock):
            targets = prefetch
        elif pressed_path is None:
            raise ValueError("Spcified prefetch but did not pass a path to pressed files")
        else:
            targets = prefetch_targets(pressed_path)
    else:
        if hmms_path is None:
            raise ValueError("Spcified prefetch but did not pass a path to the .hmm file")
        targets = pyhmmer.plan7.HMMFile(hmms_path)

    # are the sequences preloaded?
    if isinstance(seqs, str):
        seqs = pyhmmer.easel.SequenceFile(seqs, format='fasta', digital=True, alphabet=pyhmmer.easel.Alphabet.amino())
    else:
        pass

    # HMMscan execution with or without saving output to file
    if hasattr(seqs, '__len__'):
        seqs_size = len(seqs)
    else:
        seqs_size = "In file, unknown length"

    if hasattr(targets, '__len__'):
        targets_size = len(targets)
    else:
        targets_size = "In file, unknown length"

    # run hmmscan or hmmsearch
    if scan:
        logger.info(f"Running hmmscan... {seqs_size} sequences against {targets_size} HMMs, using {cpu} CPUs, additional kwargs: {kwargs}")
        all_hits = pyhmmer.hmmer.hmmscan(seqs, targets, cpus=cpu, incE=eval_con, **kwargs)
    else:
        logger.info(f"Running hmmsearch... {targets_size} HMMs against {seqs_size} seqs, using {cpu} CPUs, additional kwargs: {kwargs}")
        all_hits = pyhmmer.hmmer.hmmsearch(targets, seqs, cpus=cpu, incE=eval_con, **kwargs)
    # check if we should save the output
    if output_file is not None:
        with open(output_file, "wb") as dst:
            for i, hits in enumerate(all_hits):
                hits.write(dst, format="domains", header=i == 0)
    return all_hits


def parse_pyhmmer(all_hits, chunk_query_ids, scanned: bool = True):
    """
    Parses the TopHit pyhmmer objects, extracting query and accession IDs, and saves them to a DataFrame.

    Args:
        all_hits (list): A list of TopHit objects from pyhmmer.
        chunk_query_ids (list): A list of query IDs from the chunk.
        scanned (bool, optional): Specifies whether the sequences were scanned or searched. Defaults to True.

    Returns:
        pandas.DataFrame: A DataFrame containing the query and accession IDs.

    Notes:
        This function iterates over each protein hit in the provided list of TopHit objects and extracts the query and accession IDs.
        The resulting query and accession IDs are then saved to a DataFrame.
        Any query IDs that are missing from the parsed hits will be added to the DataFrame with a placeholder value indicating no accession information.
    """
    # initialize an empty dictionary to store the data
    parsed_hits = {}

    # iterate over each protein hit
    for top_hits in all_hits:
        for hit in top_hits:
            # check e-value
            if hit.evalue > top_hits.incE:
                continue
            # extract the query and accession IDs and decode the query ID
            if scanned:
                query_id = hit.hits.query_name.decode('utf-8')
                accession_id = hit.accession.decode('utf-8')
            else:
                query_id = hit.name.decode('utf-8')
                accession_id = hit.hits.query_accession.decode('utf-8')

            # if the query_id already exists in the dictionary, append the accession_id
            # to the existing value
            if query_id in parsed_hits:
                parsed_hits[query_id].append(accession_id)
            # otherwise, create a new key-value pair in the dictionary
            else:
                parsed_hits[query_id] = [accession_id]

    # find the query IDs that are missing from the parsed hits
    missing_query_ids = set(chunk_query_ids) - set(parsed_hits.keys())

    # add the missing query IDs with a placeholder value to indicate no
    # accession information
    for missing_query_id in missing_query_ids:
        parsed_hits[missing_query_id] = [""]

    # create the DataFrame from the dictionary
    df = pd.DataFrame(
        parsed_hits.items(),
        columns=[
            "query_id",
            "accession_id"])

    # convert list of accession IDs to string
    df["accession_id"] = df["accession_id"].apply(
        lambda x: ";".join(x) if x else "")

    return df


def local_hmmer_wrapper(chunk_index, chunked_inputs, press_path, hmm_path, out_dir, e_value: float=1e-6, prefetch=True, cpu=1, wakeup=None, scan=True, **kwargs):
    """
    A wrapping function that runs and parses pyhmmer in chunks.

    Args:
        chunk_index (int): Number of sequence chunks.
        chunked_inputs (pandas.DataFrame): DataFrame containing chunked PID inputs
        press_path (str): Path to the pressed HMMs.
        hmm_path (str): Path to the HMMs.
        out_dir (str): Path to the output directory.
        e_value (float, optional): E-value threshold. Defaults to 1e-6.
        prefetch (bool, optional): Specifies whether to prefetch the HMMs. Defaults to True.
        cpu (int, optional): Number of CPUs to use. Defaults to 1.
        wakeup (int or None, optional): Delay in seconds before starting the execution. Default is None.
        scan (bool, optional): Specifies whether to run hmmscan or hmmsearch. Defaults to True.

    Returns:
        None

    Notes:
        This function performs the following steps:
        1. Converts string sequences to pyhmmer digital blocks.
        2. Runs HMMER via pyhmmer with the provided sequences.
        3. Parses the pyhmmer output and saves it to a CSV file.

        The parsed pyhmmer output is saved in the directory specified by OUTPUT_DIR,
        with each chunk having its own separate output file named '{chunk_index}_output.csv'.

        If the wakeup parameter is specified, the function will wait for the specified
        number of seconds before starting the execution.
    """
    # we want to wait for execution to see if this worker is actually being used
    # or if it is in the process of being killed
    if wakeup is not None:
        time.sleep(wakeup)

    # convert string sequences to pyhmmer digital blocks
    sequences = save_to_digital_sequences(chunked_inputs)

    # run HMMER via pyhmmer
    if prefetch:
        hits = run_pyhmmer(
            seqs=sequences,
            pressed_path=press_path,
            prefetch=prefetch,
            cpu=cpu,
            eval_con=e_value,
            scan=scan,
            **kwargs
        )
    else:
        hits = run_pyhmmer(
            seqs=sequences,
            hmms_path=hmm_path,
            prefetch=False,
            cpu=cpu,
            eval_con=e_value,
            scan=scan,
            **kwargs
        )

    # Parse pyhmmer output and save to CSV file
    accessions_parsed = parse_pyhmmer(
        all_hits=hits, chunk_query_ids=chunked_inputs['pid'].tolist(), scanned=scan)
    accessions_parsed.to_csv(
        f'{out_dir}/{chunk_index}_output.csv',
        index=False)


#################
#### parsing ####
#################


def preprocess_accessions(meso_accession: str, thermo_accession: str):
    """
    Preprocesses meso_accession and thermo_accession by converting them to sets.

    Args:
        meso_accession (str): Meso accession string separated by ';'.
        thermo_accession (str): Thermo accession string separated by ';'.

    Returns:
        tuple: A tuple containing the preprocessed meso_accession and thermo_accession sets.
    """
    # Convert accessions to sets
    if pd.isna(meso_accession):
        meso_accession_set = set()
    else:
        meso_accession_set = set(str(meso_accession.split(';')))
    if pd.isna(thermo_accession):
        thermo_accession_set = set()
    else:
        thermo_accession_set = set(str(thermo_accession.split(';')))
    
    return meso_accession_set, thermo_accession_set


def calculate_jaccard_similarity(meso_accession_set, thermo_accession_set):
    """
    Calculates the Jaccard similarity between meso_pid and thermo_pid pairs based on their accessions.

    Jaccard similarity is defined as the size of the intersection divided by the size of the union of two sets.

    Args:
        meso_accession_set (set): Set of meso_pid accessions.
        thermo_accession_set (set): Set of thermo_pid accessions.

    Returns:
        float: Jaccard similarity between the two sets of accessions. Returns 0 if the union is empty.
    """
    # Calculate Jaccard similarity
    intersection = len(meso_accession_set.intersection(thermo_accession_set))
    union = len(meso_accession_set.union(thermo_accession_set))
    jaccard_similarity = intersection / union if union > 0 else 0

    return jaccard_similarity


def process_pairs_table(
        conn,
        dbname,
        chunk_size: int,
        output_directory,
        jaccard_threshold):
    """
    Processes the pairs table, calculates Jaccard similarity, and generates output CSV.

    Parameters:
        conn (**): Path to the database file.
        dbname (str): Name of the database.
        chunk_size (int): Size of each query chunk to fetch from the database.
        output_directory (str): Directory path to save the output CSV files.
        jaccard_threshold (float): Threshold value for Jaccard similarity.

    Returns:
        None
    """
    # Perform a join to get relevant information from the two tables
    query1 = f"""
        CREATE OR REPLACE TEMP TABLE joined_pairs AS
        SELECT p.meso_pid, p.thermo_pid, pr.accession AS meso_accession, pr2.accession AS thermo_accession
        FROM {dbname}.pairpro.final AS p
        INNER JOIN proteins_from_pairs AS pr ON (p.meso_pid = pr.pid)
        INNER JOIN proteins_from_pairs AS pr2 ON (p.thermo_pid = pr2.pid)
    """
    conn.execute(query1)

    # Define the evaluation function for the apply function
    def evaluation_function(row, jaccard_threshold):
        """
        Evaluates the Jaccard similarity between meso_pid and thermo_pid pairs based on their accessions.

        Notes:
            This function is used in the apply function to calculate the Jaccard similarity
            between meso_pid and thermo_pid pairs based on their accessions.
            There is a parsing logic for the accessions, which is described below.
            If both meso_accession and thermo_accession are nan, then the Jaccard similarity is None.
            If either meso_accession or thermo_accession is empty, then the Jaccard similarity is 0.
            If both meso_accession and thermo_accession are not empty, then the Jaccard similarity is calculated.
        """
        # Get the accessions
        meso_acc = row['meso_accession']
        thermo_acc = row['thermo_accession']
        

        # preprocessing accessions logic
        if type(meso_acc) == str:
            meso_acc_set, thermo_acc_set = preprocess_accessions(meso_acc, thermo_acc)
        elif type(meso_acc) == list:
            meso_acc_set = set(meso_acc)
            thermo_acc_set = set(thermo_acc)
        elif pd.isna(meso_acc) or pd.isna(thermo_acc):
            meso_acc_set = set()
            thermo_acc_set = set()
        else:
            raise ValueError("meso_acc must be either a string or a list")
        
        # parsing accessions logic
        if not meso_acc_set and not thermo_acc_set:
            score = None
            functional = None
        elif meso_acc_set and thermo_acc_set:
            # Preprocess the accessions
            score = calculate_jaccard_similarity(meso_acc_set, thermo_acc_set)
            functional = score > jaccard_threshold
        else:
            # Handle unmatched rows
            score = 0.0
            functional = False
            
        return {'functional': functional, 'score': score}

    # Generate output CSV file
    try:
        # Execute the query
        query = conn.execute("SELECT * FROM joined_pairs")
        data_remaining = True
        chunk_counter = 0  # Initialize the chunk counter
        while data_remaining:
            # Fetch the query result in chunks
            query_chunk = query.fetch_df_chunk(vectors_per_chunk=chunk_size)

            # Check if there is data remaining
            if len(query_chunk) == 0:
                data_remaining = False
                break

            # Calculate Jaccard similarity and determine functional status
            # using apply function
            query_chunk[['functional', 'score']] = query_chunk.apply(
                evaluation_function, axis=1, args=(jaccard_threshold,), result_type='expand')

            # Write DataFrame to CSV
            chunk_counter += 1  # Increment the chunk counter
            query_chunk.to_csv(
                f'{output_directory}{chunk_counter}_output.csv',
                index=False,
                columns=[
                    'meso_pid',
                    'thermo_pid',
                    'functional',
                    'score'])
            logger.info(f'Chunk {chunk_counter} of size {len(query_chunk)} written to csv.')

    except IOError as e:
        logger.warning(f"Error writing to CSV file: {e}")


##########################################################################################
### Variations of the above code ##########################################################
##########################################################################################


def local_hmmer_wrapper_example(chunk_index, dbpath, chunked_pid_inputs,
                                press_path, out_dir, wakeup=None):
    """
    A wrapping function that runs and parses pyhmmer in chunks.

    Args:
        chunk_index (int): Number of sequence chunks.
        dbpath (stf): Path to the database.
        chunked_pid_inputs (pandas.DataFrame): DataFrame containing chunked PID inputs.
        press_path (str): Path to the pressed HMM database.
        out_path (str): Path to directory where output will be saved.
        wakeup (int or None, optional): Delay in seconds before starting the execution. Default is None.

    Returns:
        None

    Notes:
        This function performs the following steps:
        1. Queries the database to get sequences only from chunked_pid_inputs.
        2. Converts the query result to a DataFrame.
        3. Converts string sequences to pyhmmer digital blocks.
        4. Runs HMMER via pyhmmer with the provided sequences.
        5. Parses the pyhmmer output and saves it to a CSV file.

        The parsed pyhmmer output is saved in the directory specified by OUTPUT_DIR,
        with each chunk having its own separate output file named '{chunk_index}_output.csv'.

        If the wakeup parameter is specified, the function will wait for the specified
        number of seconds before starting the execution.
    """
    # we want to wait for execution to see if this worker is actually being used
    # or if it is in the process of being killed
    if wakeup is not None:
        time.sleep(wakeup)

    # query the database to get sequences only from chunked_pid_inputs
    conn = ddb.connect(dbpath, read_only=True)

    # get the unique pids from the chunked_pid_inputs
    pids = set(chunked_pid_inputs["pid"])

    # Only extract protein_seqs from the list of PID inputs
    placeholders = ', '.join(['?'] * len(pids))
    query = f"SELECT pid, protein_seq FROM pairpro.proteins WHERE pid IN ({placeholders})"
    query_db = conn.execute(query, list(pids)).fetchall()

    # close db connection
    conn.close()

    # convert the query db to a dataframe
    result_df = pd.DataFrame(query_db, columns=['pid', 'protein_seq'])

    # convert string sequences to pyhmmer digital blocks
    sequences = save_to_digital_sequences(result_df)

    # run HMMER via pyhmmer
    hits = run_pyhmmer(
        seqs=sequences,
        hmms_path=press_path,
        press_path=press_path,
        prefetch=True,
        cpu=1,
        eval_con=1e-12)

    # get the query IDs from the chunked_pid_inputs
    chunk_query_ids = chunked_pid_inputs["pid"].tolist()

    # Parse pyhmmer output and save to CSV file
    accessions_parsed = parse_pyhmmer(
        all_hits=hits, chunk_query_ids=chunk_query_ids)
    accessions_parsed.to_csv(
        f'{out_dir}/{chunk_index}_output.csv',
        index=False)

def process_pairs_table_ana(
        conn,
        dbname,
        chunk_size: int,
        output_directory,
        jaccard_threshold):
    """
    Processes the pairs table, calculates Jaccard similarity, and generates output CSV.

    Parameters:
        conn (**): Path to the database file.
        dbname (str): Name of the database.
        chunk_size (int): Size of each query chunk to fetch from the database.
        output_directory (str): Directory path to save the output CSV files.
        jaccard_threshold (float): Threshold value for Jaccard similarity.

    Returns:
        None
    """
    # Perform a join to get relevant information from the two tables
    query1 = f"""
        CREATE OR REPLACE TEMP TABLE joined_pairs AS
        SELECT p.meso_pid, p.thermo_pid, pr.accession AS meso_accession, pr2.accession AS thermo_accession
        FROM {dbname}.pairpro.final AS p
        INNER JOIN proteins_from_pairs4 AS pr ON (p.meso_pid = pr.pid)
        INNER JOIN proteins_from_pairs4 AS pr2 ON (p.thermo_pid = pr2.pid)
    """
    conn.execute(query1)

    # Define the evaluation function for the apply function
    def evaluation_function(row, jaccard_threshold):
        """
        Evaluates the Jaccard similarity between meso_pid and thermo_pid pairs based on their accessions.

        Notes:
            This function is used in the apply function to calculate the Jaccard similarity
            between meso_pid and thermo_pid pairs based on their accessions.
            There is a parsing logic for the accessions, which is described below.
            If both meso_accession and thermo_accession are nan, then the Jaccard similarity is None.
            If either meso_accession or thermo_accession is empty, then the Jaccard similarity is 0.
            If both meso_accession and thermo_accession are not empty, then the Jaccard similarity is calculated.
        """
        # Get the accessions
        meso_acc = row['meso_accession']
        thermo_acc = row['thermo_accession']
        
        # parsing accessions logic
        if meso_acc == 'nan' and thermo_acc == 'nan':
            score = None
            functional = None
        elif meso_acc and thermo_acc:
            # Preprocess the accessions
            meso_acc_set, thermo_acc_set = preprocess_accessions(
                meso_acc, thermo_acc)
            score = calculate_jaccard_similarity(meso_acc_set, thermo_acc_set)
            functional = score > jaccard_threshold
        else:
            # Handle unmatched rows
            score = 0.0
            functional = False

        return {'functional': functional, 'score': score}

    # Generate output CSV file
    try:
        # Execute the query
        query = conn.execute("SELECT * FROM joined_pairs")
        data_remaining = True
        chunk_counter = 0  # Initialize the chunk counter
        while data_remaining:
            # Fetch the query result in chunks
            query_chunk = query.fetch_df_chunk(vectors_per_chunk=chunk_size)

            # Check if there is data remaining
            if len(query_chunk) == 0:
                data_remaining = False
                break

            # Calculate Jaccard similarity and determine functional status
            # using apply function
            query_chunk[['functional', 'score']] = query_chunk.apply(
                evaluation_function, axis=1, args=(jaccard_threshold,), result_type='expand')

            # Write DataFrame to CSV
            chunk_counter += 1  # Increment the chunk counter
            query_chunk.to_csv(
                f'{output_directory}{chunk_counter}_output.csv',
                index=False,
                columns=[
                    'meso_pid',
                    'thermo_pid',
                    'functional',
                    'score'])
            logger.info(f'Chunk {chunk_counter} of size {len(query_chunk)} written to csv.')

    except IOError as e:
        logger.warning(f"Error writing to CSV file: {e}")