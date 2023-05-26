"""
The following is importable code for running HMMER 
either locally via pyhmmer or via Interpro's API

TODO:
- Amin update documentation
"""
# system dependecies
import asyncio
from concurrent.futures import ProcessPoolExecutor
import json
import logger
import nest_asyncio
import os
import requests
import time
import urllib.parse

# library dependencies
import duckdb as ddb
import httpx
from joblib import Parallel, delayed
import pandas as pd
import pyhmmer

# local dependencies
import FAFSA.utils

####### API HMMER
async def send_request(semaphore, sequence, client):
    """
    Sends a POST request to the HMMER API with a protein sequence.
    -------------
    Parameters:
    -------------
    semaphore: asyncio.Semaphore
        A semaphore to limit concurrent requests.
    sequence: str
        The protein sequence to be sent in the request.
    client: httpx.AsyncClient
        An HTTP client for sending the request.
    -------------
    Returns:
    -------------
    response: httpx.Response
        The response received from the HMMER API.
    """

    url = 'https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan'
    headers = {'Content-Type': 'application/x-www-form-urlencoded',
               'Accept': 'application/json'}
    data = {'hmmdb': 'pfam', 'seq': f'>seq\n{sequence}'}
    data = urllib.parse.urlencode(data).encode('ascii')

    async with semaphore:
        response = await client.post(url, headers=headers, data=data, follow_redirects=False, timeout=15000)

    return response


async def process_response(semaphore, sequence, response, client, pid, max_retries=3):
    """
    Processes the response received from the HMMER API.
    -------------
    Parameters:
    -------------
    semaphore: asyncio.Semaphore
        A semaphore to limit concurrent requests.
    sequence: str
        The protein sequence associated with the response.
    response: httpx.Response
        The response received from the HMMER API.
    client: httpx.AsyncClient
        An HTTP client for sending subsequent requests.
    pid: int
        The protein ID associated with the sequence.
    max_retries: int, optional
        The maximum number of retries for failed requests (default is 3).
    -------------
    Returns:
    -------------
    dfff: pd.DataFrame or None
        A DataFrame containing the search results for the protein sequence, or None if an error occurred.
    """

    redirect_url = response.headers.get('Location')

    if redirect_url is None:
        print("Error: No redirect URL found in response.")
    else:
        headers = {'Accept': 'application/json'}
        async with semaphore:
            for attempt in range(max_retries):
                try:
                    response2 = await client.get(redirect_url, headers=headers, timeout=15000)
                    break
                except httpx.ReadTimeout:
                    if attempt < max_retries - 1:
                        # Exponential backoff
                        await asyncio.sleep(5 ** attempt)
                    else:
                        raise
        try:
            results = response2.json()
            hits = results['results']['hits']
        except KeyError:
            logger.info(
                f"Error: 'results' key not found in response for sequence {sequence}.")
            return None
        except json.JSONDecodeError:
            logger.info(
                f"Error: JSONDecodeError for sequence {sequence}. Response text: {response2.text}")
            return None

        if hits:
            loop = asyncio.get_event_loop()
            dfff = await loop.run_in_executor(None, pd.json_normalize, hits, 'domains', ['acc', 'name', 'score', 'evalue', 'pvalue', 'desc'])
            dfff.insert(0, 'sequence', sequence)
            # Add new column here
            dfff.insert(0, 'pid', pid)
            dfff = dfff.set_index('pid')  # Set new column as index
            return dfff
        else:
            return None


async def hmmerscanner(df: pd.DataFrame, k: int, max_concurrent_requests: int, output_path: str):
    """
    Runs the HMMER scanner for protein sequences.
    -------------
    Parameters:
    -------------
    df: pd.DataFrame
        A DataFrame that contains protein sequences.
    k: int
        The number of protein sequences to search.
    max_concurrent_requests: int
        The maximum number of concurrent requests to the HMMER API.
    output_path: str
        The output directory of where the data will be stored.
    -------------
    Returns:
    -------------
    results_df: pd.DataFrame
        A DataFrame containing the search results for all protein sequences.
    """

    if k > 1000:
        print("Use local function for the number of sequences more than 1000.")
        return pd.DataFrame()

    sequences = df['protein_seq'][:k]
    # Get corresponding prot_pair_index values
    indices = df['pid'][:k]
    tasks = []
    semaphore = asyncio.Semaphore(max_concurrent_requests)

    # Use a process pool to parallelize JSON processing and DataFrame creation
    with ProcessPoolExecutor() as executor:
        loop = asyncio.get_event_loop()
        async with httpx.AsyncClient() as client:
            for seq, idx in zip(sequences, indices):  # Include the index here
                task = asyncio.create_task(
                    send_request(semaphore, seq, client))
                tasks.append(task)

            responses = await asyncio.gather(*tasks)

            tasks = []
            for (seq, idx), response in zip(zip(sequences, indices), responses):  # Include the index here
                task = asyncio.create_task(process_response(
                    semaphore, seq, response, client, idx))  # idx is the prot
                tasks.append(task)

            results = await asyncio.gather(*tasks)
    common_columns = set.intersection(
        *(set(df.columns) for df in results if df is not None))
    results_df = pd.concat(
        [result[list(common_columns)] for result in results if result is not None])
    # write result to csv
    results_df.to_csv(f'{output_path}API_output.csv')
    return results_df


def run_hmmerscanner(df: pd.DataFrame, k: int, max_concurrent_requests: int):
    """
    This function runs the hmmerscanner function within an event loop and returns the search results as
    a DataFrame.
    -------------
    Parameters:
    -------------
    df: pandas.core.DataFrame
    A DataFrame that contains protein sequences.
    k: int
    The number of protein sequences to search.
    max_concurrent_requests: int
    The maximum number of concurrent requests to the HMMER API.
    -------------
    Returns:
    -------------
    
    results_df: pandas.core.DataFrame
        A DataFrame containing the search results for all protein sequences.
    """

    # Set up the event loop and call the hmmerscanner function
    nest_asyncio.apply()
    return asyncio.run(hmmerscanner(df, k, max_concurrent_requests))


####### Local HMMER

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
        Path to the pressed HMM database.

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


def save_to_digital_sequences(dataframe: pd.DataFrame):
    """
    Save protein sequences from a DataFrame to a digital sequence block.

    Parameters
    ----------
    dataframe : pd.DataFrame
        DataFrame containing PIDs (Protein IDs) and sequences.

    Returns
    -------
    DigitalSequenceBlock
        A digital sequence block containing the converted sequences.
    """
    # Create empty list
    seqlist = []

    # Establish pyhmmer alphabet
    amino_acids = pyhmmer.easel.Alphabet.amino()

    # Convert proteins in dataframe to suitable format
    for _, row in dataframe.iterrows():
        pid = bytes(row['pid'], encoding='utf-8')
        seq_str = row['protein_seq']
        sequences = pyhmmer.easel.TextSequence(name=pid, sequence= seq_str)
        sequences = sequences.digitize(amino_acids)
        seqlist.append(sequences)
    
    # Convert so SequenceBlocks
    seqblock = pyhmmer.easel.DigitalSequenceBlock(amino_acids, seqlist)

    return seqblock


def run_pyhmmer(
        seqs: pyhmmer.easel.DigitalSequenceBlock,
        hmms_path: str,
        press_path: str,
        prefetch: bool = False,
        output_file: str = None,
        cpu: int = 4,
        eval_con: float = 1e-10):
    """
    Run HMMER's hmmscan program on a set of input sequences using with HMMs from a database.
    Parameters
    ----------
    seqs : pyhmmer.easel.DigitalSequenceBlock
        Path to the input sequence file.

    hmms_path : str
        Path to the HMM database.

    press_path : str
        Path to the pressed HMM database

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
    # ensure output_file has .domtblout extension
    if output_file is not None and not output_file.endswith('.domtblout'):
        output_file = f"{os.path.splitext(output_file)[0]}.domtblout"


    # HMM profile modes
    if prefetch:
        targets = prefetch_targets(hmms_path)
    else:
        targets = pyhmmer.plan7.HMMFile(press_path)

    # HMMscan execution with or without saving output to file
    all_hits = list(pyhmmer.hmmer.hmmscan(seqs, targets, cpus=cpu, E=eval_con))
    # check if we should save the output
    if output_file is not None:
        with open(output_file, "wb") as dst:
            for i, hits in enumerate(all_hits):
                hits.write(dst, format="domains", header=i == 0)
    return all_hits


def parse_pyhmmer(all_hits, chunk_query_ids):
    """
    Parses the TopHit pyhmmer object getting the query and accession IDs and saves to a DataFrame

    Parameters
    ----------
    all_hits : list
        A list of TopHit objects from pyhmmer.
    chunk_query_ids : list
        A list of query IDs from the chunk.

    Returns
    -------
    pandas.DataFrame
        A dataframe containing the query and accession IDs.
    """
    # initialize an empty dictionary to store the data
    parsed_hits = {}

    # iterate over each protein hit
    for top_hits in all_hits:
        for hit in top_hits:
            # extract the query and accession IDs and decode the query ID
            query_id = hit.hits.query_name.decode('utf-8')
            accession_id = hit.accession.decode('utf-8')

            # if the query_id already exists in the dictionary, append the accession_id
            # to the existing value
            if query_id in parsed_hits:
                parsed_hits[query_id].append(accession_id)
            # otherwise, create a new key-value pair in the dictionary
            else:
                parsed_hits[query_id] = [accession_id]

    # find the query IDs that are missing from the parsed hits
    missing_query_ids = set(chunk_query_ids) - set(parsed_hits.keys())

    # add the missing query IDs with a placeholder value to indicate no accession information
    for missing_query_id in missing_query_ids:
        parsed_hits[missing_query_id] = ['No Accession Information']

    # create the DataFrame from the dictionary and convert list of accession IDs to string
    df = pd.DataFrame(parsed_hits.items(), columns=["query_id", "accession_id"])
    df["accession_id"] = df["accession_id"].apply(lambda x: ';'.join(x))

    return df


def local_hmmer_wrapper(chunk_index, dbpath, chunked_pid_inputs, 
                        press_path, out_dir,  wakeup=None):
    """
    A wrapping function that runs and parses pyhmmer in chunks.

    Parameters
    ----------
    chunk_index : int
        Number of sequence chunks.

    dbpath : str
        Path to the database.

    chunked_pid_inputs : pandas.DataFrame
        Dataframe containing chunked PID inputs.

    press_path : str
        Path to the pressed HMM database

    out_path : str
        Path to directory where output will be saved

    wakeup : int or None, optional
        Delay in seconds before starting the execution, by default None.
    Returns
    -------
    None

    Notes
    -----
    This function performs the following steps:
    1. Queries the database to get sequences only from chunked_pid_inputs.
    2. Converts the query result to a dataframe.
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
    query = f"SELECT pid, protein_seq FROM proteins WHERE pid IN ({placeholders})"
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
        prefetch=True,
        cpu=1,
        eval_con=1e-5)
    
    # get the query IDs from the chunked_pid_inputs
    chunk_query_ids = chunked_pid_inputs["pid"].tolist()

    # Parse pyhmmer output and save to CSV file
    accessions_parsed = parse_pyhmmer(all_hits=hits, chunk_query_ids=chunk_query_ids)
    accessions_parsed.to_csv(
        f'{out_dir}/{chunk_index}_output.csv',
        index=False)