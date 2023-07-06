"""
Adapatations of different functions for user input scenarios.
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


# API HMMER


async def send_request(semaphore, sequence, client):
    """
    Asynchronously sends a POST request to the HMMER API, submitting a protein sequence for analysis.

    Args:
        semaphore (asyncio.Semaphore): An object that controls concurrent request submission, helping to avoid server overload.
        sequence (str): The protein sequence that is to be analyzed and included in the body of the POST request.
        client (httpx.AsyncClient): An HTTP client for sending the request.

    Returns:
        httpx.Response: Response received from the HMMER API.

    Raises:
        httpx.HTTPStatusError: If the HTTP request returned a status code that denotes an error.
        httpx.TimeoutException: If the request times out.
    """

    url = 'https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan'
    headers = {'Content-Type': 'application/x-www-form-urlencoded',
               'Accept': 'application/json'}
    data = {'hmmdb': 'pfam', 'seq': f'>seq\n{sequence}'}
    data = urllib.parse.urlencode(data).encode('ascii')

    async with semaphore:
        response = await client.post(url, headers=headers, data=data, follow_redirects=False, timeout=15000)

    return response


async def process_response(semaphore, sequence, response, client, pair_id, max_retries=3):
    """
    Processes the response received from the HMMER API, including retrying requests that have failed.

    Args:
        semaphore (asyncio.Semaphore): An object that controls concurrent request submission, helping to avoid server overload.
        sequence (str): The protein sequence associated with the response.
        response (httpx.Response): The response received from the HMMER API.
        client (httpx.AsyncClient): An HTTP client for sending subsequent requests.
        pair_id (int): The protein ID associated with the sequence.
        max_retries (int, optional): The maximum number of retries for failed requests. Defaults to 3.

    Returns:
        pd.DataFrame or None: A DataFrame containing the search results for the protein sequence,
        or None if an error occurred.

    Raises:
        KeyError: If expected key is not found in the response.
        json.JSONDecodeError: If JSON decoding fails.
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
            dfff.insert(0, 'pair_id', pair_id)
            dfff = dfff.set_index('pair_id')  # Set new column as index
            return dfff
        else:
            return None


async def hmmerscanner(df: pd.DataFrame, which: str, k: int, max_concurrent_requests: int, output_path: str):
    """
    Scans multiple protein sequences using the HMMER API, asynchronously submitting and processing each request.

    Args:
        df (pd.DataFrame): A DataFrame containing protein sequences.
        which (str): The column name of the protein sequences.
        k (int): The number of protein sequences to search.
        max_concurrent_requests (int): The maximum number of concurrent requests to the HMMER API.
        output_path (str): The output directory where the data will be stored.

    Returns:
        pd.DataFrame: A DataFrame containing the search results for all protein sequences.

    Raises:
        ValueError: If the number of sequences exceeds the limit of 1000.
    """

    if k > 1000:
        print("Use local function for the number of sequences more than 1000.")
        return pd.DataFrame()

    sequences = df[which][:k]
    # Get corresponding prot_pair_index values
    indices = df['pair_id'][:k]
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
            for (seq, idx), response in zip(
                    zip(sequences, indices), responses):  # Include the index here
                task = asyncio.create_task(process_response(
                    semaphore, seq, response, client, idx))  # idx is the prot
                tasks.append(task)

            results = await asyncio.gather(*tasks)
    common_columns = list(set.intersection(
        *(set(df.columns) for df in results if df is not None)))
    results_df = pd.concat([result[list(common_columns)]
                            for result in results if result is not None])
    # write result to csv
    results_df.to_csv(f'{output_path}{which}_API_output.csv')
    return results_df


def run_hmmerscanner(
        df: pd.DataFrame,
        which: str,
        k: int,
        max_concurrent_requests: int,
        output_path: str):
    """
    Runs the asynchronous HMMER scanning operation in a new event loop.

    Args:
        df (pd.DataFrame): A DataFrame containing protein sequences.
        which (str): The column name of the protein sequences.
        k (int): The number of protein sequences to search.
        max_concurrent_requests (int): The maximum number of concurrent requests to the HMMER API.
        output_path (str): The output directory where the data will be stored. (like '/Users/amin/ValidProt/data/')

    Returns:
        pd.DataFrame: A DataFrame containing the search results for all protein sequences.

    Raises:
        nest_asyncio.NestingError: If the event loop is already running.
        Any exceptions raised by hmmerscanner function.
    """

    # Set up the event loop and call the hmmerscanner function
    nest_asyncio.apply()
    return asyncio.run(hmmerscanner(df, which, k, max_concurrent_requests, output_path))


# user input


def save_to_digital_sequences_user_query(dataframe: pd.DataFrame):
    """
    Save protein sequences from a DataFrame to a digital sequence block.

    Args:
        dataframe (pd.DataFrame): DataFrame containing pair_id (Protein pair IDs) and sequences.

    Returns:
        pyhmmer.easel.DigitalSequenceBlock: A digital sequence block containing the converted sequences.
    """
    # Create empty list
    query_seqlist = []

    # Establish pyhmmer alphabet
    amino_acids = pyhmmer.easel.Alphabet.amino()

    # Convert proteins in dataframe to suitable format
    for _, row in dataframe.iterrows():
        pair_id = bytes(str(row['pair_id']), encoding='utf-8')
        seq_str = row['query']
        sequences = pyhmmer.easel.TextSequence(name=pair_id, sequence=seq_str)
        sequences = sequences.digitize(amino_acids)
        query_seqlist.append(sequences)

    # Convert so SequenceBlocks
    query_seqblock = pyhmmer.easel.DigitalSequenceBlock(
        amino_acids, query_seqlist)

    return query_seqblock


def save_to_digital_sequences_user_subject(dataframe: pd.DataFrame):
    """
    Save protein sequences from a DataFrame to a digital sequence block.

    Args:
        dataframe (pd.DataFrame): DataFrame containing pair_id (Protein pair IDs) and sequences.

    Returns:
        pyhmmer.easel.DigitalSequenceBlock: A digital sequence block containing the converted sequences.
    """
    # Create empty list
    subject_seqlist = []

    # Establish pyhmmer alphabet
    amino_acids = pyhmmer.easel.Alphabet.amino()

    # Convert proteins in dataframe to suitable format
    for _, row in dataframe.iterrows():
        pair_id = bytes(str(row['pair_id']), encoding='utf-8')
        seq_str = row['subject']
        sequences = pyhmmer.easel.TextSequence(name=pair_id, sequence=seq_str)
        sequences = sequences.digitize(amino_acids)
        subject_seqlist.append(sequences)

    # Convert so SequenceBlocks
    subject_seqblock = pyhmmer.easel.DigitalSequenceBlock(
        amino_acids, subject_seqlist)

    return subject_seqblock


def parse_pyhmmer_user(all_hits, chunk_pair_ids):
    """
    Parses the TopHit pyhmmer objects, extracting query and accession IDs, and saves them to a DataFrame.

    Args:
        all_hits (list): A list of TopHit objects from pyhmmer.
        chunk_pair_ids (list): A list of query IDs from the chunk.

    Returns:
        pandas.DataFrame: A DataFrame containing the pair and accession IDs.

    Notes:
        This function iterates over each protein hit in the provided list of TopHit objects and extracts the query and accession IDs.
        The resulting query and accession IDs are then saved to a DataFrame.
        Any pair IDs that are missing from the parsed hits will be added to the DataFrame with a placeholder value indicating no accession information.
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
    missing_query_ids = set(chunk_pair_ids) - set(parsed_hits.keys())

    # add the missing query IDs with a placeholder value to indicate no
    # accession information
    for missing_query_id in missing_query_ids:
        parsed_hits[missing_query_id] = [""]

    # create the DataFrame from the dictionary
    df = pd.DataFrame(parsed_hits.items(), columns=["pair_id", "accession_id"])

    # convert list of accession IDs to string
    df["accession_id"] = df["accession_id"].apply(
        lambda x: ";".join(x) if x else "")

    return df


def user_local_hmmer_wrapper_query(
        chunk_index,
        press_path,
        sequences,
        out_dir):
    """
    TODO
    """
    # convert string sequences to pyhmmer digital blocks
    query_seqblock = save_to_digital_sequences_user_query(sequences)

    # run HMMER via pyhmmer
    query_hits = run_pyhmmer(
        seqs=query_seqblock,
        hmms_path=press_path,
        press_path=press_path,
        prefetch=True,
        cpu=1,
        eval_con=1e-5)

    # get the query IDs from the chunked_pid_inputs
    missing_pair_ids = sequences["pair_id"].tolist()

    # Parse pyhmmer output and save to CSV file
    accessions_parsed_query = parse_pyhmmer_user(
        all_hits=query_hits, chunk_pair_ids=missing_pair_ids)
    accessions_parsed_query.to_csv(
        f"{out_dir}query_result_{chunk_index}.csv",
        index=False)


def user_local_hmmer_wrapper_subject(
        chunk_index,
        press_path,
        sequences,
        out_dir):
    """
    TODO
    """
    # convert string sequences to pyhmmer digital blocks
    subject_seqblock = save_to_digital_sequences_user_subject(sequences)

    # run HMMER via pyhmmer
    subject_hits = run_pyhmmer(
        seqs=subject_seqblock,
        hmms_path=press_path,
        press_path=press_path,
        prefetch=True,
        cpu=1,
        eval_con=1e-5)

    # get the query IDs from the chunked_pid_inputs
    missing_pair_ids = sequences["pair_id"].tolist()

    accessions_parsed_subject = parse_pyhmmer_user(
        all_hits=subject_hits, chunk_pair_ids=missing_pair_ids)
    accessions_parsed_subject.to_csv(
        f"{out_dir}subject_result_{chunk_index}.csv",
        index=False)


#### API parsing #####
def parse_function_csv_API(file_path: str) -> Dict[str, List[str]]:
    """
    Parses the CSV file with protein IDs and their corresponding accession IDs
    and returns a dictionary with protein IDs as keys and accession IDs as values.
    """
    # Create a dictionary to store csv results
    protein_dict = defaultdict(list)

    with open(file_path, 'r') as csvfile:
        # read csv
        reader = csv.reader(csvfile)
        # get the header row
        header = next(reader)
        # find the index of the pair_id and accession columns
        pair_id_idx = header.index('pair_id')
        accession_idx = header.index('acc')

        for row in reader:
            pair_id = row[pair_id_idx]
            accessions = row[accession_idx].split(';')
            protein_dict[pair_id].extend(accessions)

    return protein_dict


def find_jaccard_similarity_API(set1: set, set2: set) -> float:
    """
    Calculates the Jaccard similarity score between two sets.
    """
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    if union == 0:
        return 0.0
    else:
        return intersection / union


def calculate_similarity_API(
        file1: str, file2: str, threshold: float) -> Dict[str, Tuple[str, float]]:
    """
    Calculates the Jaccard similarity score between each protein in file1 and file2,
    and returns a dictionary with query IDs as keys and a tuple indicating whether
    the score threshold was met and the Jaccard similarity score.
    """
    # Read the CSV files and create dictionaries with query IDs and accession
    # IDs
    dict1 = parse_function_csv_API(file1)
    dict2 = parse_function_csv_API(file2)

    # Create a dictionary to store the Jaccard similarity scores
    scores = defaultdict(float)

    # Calculate the Jaccard similarity score between each protein in file1 and
    # file2
    for query1, accs1 in dict1.items():
        for query2, accs2 in dict2.items():
            if query1 == query2:
                score = find_jaccard_similarity_API(set(accs1), set(accs2))
                scores[(query1, query2)] = score

    # Create a dictionary to store the functional tuple values
    functional = {}

    # Set the functional tuple value based on the Jaccard similarity score
    # threshold
    for (query1, query2), score in scores.items():
        if score >= threshold:
            functional[(query1, query2)] = (True, score)
        else:
            functional[(query1, query2)] = (False, score)

    return functional


def write_function_output_API(
        output_dict: Dict[str, Tuple[str, float]], output_file: str):
    """
    Writes a dictionary of protein pair IDs and functional tuple values to a CSV file.

    Args:
    output_dict : Dict[str, Tuple[str, float]]
        A dictionary of protein pair IDs and functional tuple values
    output_file : str
        File path to write the output CSV file
    """
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['file1', 'file2', 'functional', 'jaccard Score']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for pair_id, (functional, score) in output_dict.items():
            writer.writerow({
                'file1': pair_id[0],
                'file2': pair_id[1],
                'functional': functional,
                'jaccard Score': score
            })


def get_file_pairs_API(directory_path):
    """
    A quick silly function to get pairs
    """
    subject_files = glob.glob(f"{directory_path}/subject_API_output.csv")
    query_files = glob.glob(f"{directory_path}/query_API_output.csv")
    print(subject_files)
    subject_files.sort()
    query_files.sort()
    file_pairs = []
    for subject_file, query_file in zip(subject_files, query_files):
        file_pairs.append((subject_file, query_file))
    return file_pairs


def parse_function_csv_user(file_path: str) -> Dict[str, List[str]]:
    """
    Parses the CSV file with protein IDs and their corresponding accession IDs
    and returns a dictionary with protein IDs as keys and accession IDs as values.
    """
    # Create a dictionary to store csv results
    protein_dict = {}

    with open(file_path, 'r') as csvfile:
        # read csv
        reader = csv.reader(csvfile)
        # skip header
        next(reader)
        for row in reader:
            query_id = row[0]
            accessions = row[1].split(';')
            protein_dict[query_id] = accessions
    return protein_dict


def calculate_similarity_user(
        file1: str, file2: str, threshold: float) -> Dict[str, Tuple[str, float]]:
    """
    Calculates the Jaccard similarity score between each protein in file1 and file2,
    and returns a dictionary with query IDs as keys and a tuple indicating whether
    the score threshold was met and the Jaccard similarity score.
    """
    # Read the CSV files and create dictionaries with query IDs and accession
    # IDs
    dict1 = parse_function_csv_user(file1)
    dict2 = parse_function_csv_user(file2)

    # Create a dictionary to store the Jaccard similarity scores
    scores = defaultdict(float)

    # Calculate the Jaccard similarity score between each protein in file1 and
    # file2
    for query1, accs1 in dict1.items():
        for query2, accs2 in dict2.items():
            if query1 == query2:
                score = find_jaccard_similarity_API(set(accs1), set(accs2))
                scores[(query1, query2)] = score

    # Create a dictionary to store the functional tuple values
    functional = {}

    # Set the functional tuple value based on the Jaccard similarity score
    # threshold
    for (query1, query2), score in scores.items():
        if score >= threshold:
            functional[(query1, query2)] = (True, score)
        else:
            functional[(query1, query2)] = (False, score)

    return functional


def get_file_pairs_user(directory_path):
    """
    A quick silly function to get pairs
    """
    query_files = glob.glob(f"{directory_path}/query_result_*.csv")
    subject_files = glob.glob(f"{directory_path}/subject_result_*.csv")
    print(query_files)
    query_files.sort()
    subject_files.sort()
    file_pairs = []
    for query_file, subject_file in zip(query_files, subject_files):
        query_chunk_index = int(query_file.split("_")[-1].split(".")[0])
        subject_chunk_index = int(subject_file.split("_")[-1].split(".")[0])
        if query_chunk_index == subject_chunk_index:
            file_pairs.append((query_file, subject_file))
    return file_pairs
