"""
This is script is mainly to give Logan a function label so he can work.
It will be reworked after the learn2thermDB script.
"""

# system dependecies
from collections import defaultdict
import csv
import glob
import os
from pathlib import Path
import sys
from typing import Dict, List, Tuple

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


def parse_function_csv(file_path: str) -> Dict[str, List[str]]:
    """
    Reads a CSV file and returns a dictionary with query IDs as keys and a list
    of accession IDs as values.

    Parameters
    ----------
    file_path : str
        File path to the HMMER csv output

    Returns
    -------
    Dict[str, List[str]]
        a dictionary consisting of the query ID and accession IDs as strings
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


def find_jaccard_similarity(set1, set2):
    """
    Calculates the Jaccard similarity score between two sets.
    """
    intersection = set1.intersection(set2)
    union = set1.union(set2)
    return len(intersection) / len(union)


def calculate_similarity(file1: str, file2: str, threshold: float) -> Dict[str, Tuple[str, float]]:
    """
    Calculates the Jaccard similarity score between each protein in file1 and file2,
    and returns a dictionary with query IDs as keys and a tuple indicating whether
    the score threshold was met and the Jaccard similarity score.
    """
    # Read the CSV files and create dictionaries with query IDs and accession IDs
    dict1 = parse_function_csv(file1)
    dict2 = parse_function_csv(file2)
    
    # Create a dictionary to store the Jaccard similarity scores
    scores = defaultdict(float)
    
    # Calculate the Jaccard similarity score between each protein in file1 and file2
    for query1, accs1 in dict1.items():
        for query2, accs2 in dict2.items():
            score = find_jaccard_similarity(set(accs1), set(accs2))
            scores[(query1, query2)] = score
    
    # Create a dictionary to store the functional tuple values
    functional = {}
    
    # Set the functional tuple value based on the Jaccard similarity score threshold
    for (query1, query2), score in scores.items():
        if score >= threshold:
            functional[(query1, query2)] = ('Yes', score)
        else:
            functional[(query1, query2)] = ('No', score)
    
    return functional



def write_function_output(output_dict: Dict[str, Tuple[str, float]], output_file: str):
    """
    Writes a dictionary of protein query IDs and functional tuple values to a CSV file.

    Parameters
    ----------
    output_dict : Dict[str, Tuple[str, float]]
        A dictionary of protein query IDs and functional tuple values
    output_file : str
        File path to write the output CSV file
    """
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['File1', 'File2', 'Functional?', 'Jaccard Score']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        for query, (functional, score) in output_dict.items():
            writer.writerow({
                'File1': query[0],
                'File2': query[1],
                'Functional?': functional,
                'Jaccard Score': score
            })


if __name__ == '__main__':
    # Logging
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    fh = logging.FileHandler('function_labelling.log', mode='w')
    formatter = logging.Formatter(
        '%(filename)-12s %(asctime)s;%(funcName)-12s: %(levelname)-8s %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    logger.info("TEST LOG")

    # set path directory
    dir_path = "/Users/humoodalanzi/ValidProt/scripts/results"

    def get_file_pairs(directory_path):
        """
        A quick silly function to get pairs
        """
        meso_files = glob.glob(f"{directory_path}/meso_result_*.csv")
        thermo_files = glob.glob(f"{directory_path}/thermo_result_*.csv")
        print(meso_files)
        meso_files.sort()
        thermo_files.sort()
        file_pairs = []
        for meso_file, thermo_file in zip(meso_files, thermo_files):
            meso_chunk_index = int(meso_file.split("_")[-1].split(".")[0])
            thermo_chunk_index = int(thermo_file.split("_")[-1].split(".")[0])
            if meso_chunk_index == thermo_chunk_index:
                file_pairs.append((meso_file, thermo_file))
        return file_pairs

    
    # set threshold
    threshold = 0.33

    logger.info("Parameters set")

    logger.info("Calculating similarity")
    # Get file pairs and calculate similarity for each pair
    file_pairs = get_file_pairs(dir_path)
    print(file_pairs)
    print(dir_path)
    results = {}
    for file1, file2 in file_pairs:
        logger.info(f"Processing {file1} and {file2}")
        output_file = f"functional_{os.path.basename(file1)}"
        similarity_scores = calculate_similarity(file1, file2, threshold)
        write_function_output(similarity_scores, output_file)
        results[(file1, file2)] = output_file

    logger.info("All pairs processed")
    logger.info("Results:")
    for pair, output_file in results.items():
        logger.info(f"{pair[0]} and {pair[1]}: {output_file}")