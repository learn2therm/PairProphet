'''
The package you need to run this script are the following:
- numpy
- biopython
'''

# system dependecies
import subprocess
from pathlib import Path


# library dependencies
from collections import defaultdict
import pandas as pd



## biopython
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SearchIO


# local dependencies/utils

## Paths
PFAM_PATH = Path("/Users/humoodalanzi/pfam/Pfam-A.hmm")
ID_DB_PATH = Path("/Users/humoodalanzi/pfam/proteins_id.zip")



def hmmer_wrapper(seq: pd.core.series.Series, input_filename: str, pfam_path: str, input_filename_with_ext: str, output_filename_with_ext: str, cpu: int = 4):
    """
    Executes the HMMER against pfam and parsing the results

    Parameters:
    ------------
    seq : str
        a pandas series of of amino acids (has to be processed in a certain way)
    input_filename : str
        A file name for the input of the transformed seq to FASTA
    pfam_path : str
        path of the HMMER/pfam db
    input_filename_with_ext : str
        A file name for the input FASTA has to include the ext. FASTA
    output_filename_with_ext : str
        output file name perferred extension is domtblout
    cpu : 4
        number of cpus for i/o
    
    Returns:
    ------------
    Input fasta file
    Output domtblout file
        
    Raises
    -------
    TypeError
    """
    
    # Create a list of SeqRecord objects
    # Write the list of SeqRecord objects to a FASTA file
    def read_seq(lists: pd.core.series.Series, inputname: str = "input"):
        """
        Returns a list of SeqRecord objects and creates a corresponding input Fasta of them

        Parameters:
        ------------
        list : pandas.core.series.Series
            a list of string amino acid sequences
        input name : str, default = 'input'
            a name for the input fasta file

        
        Returns:
        ------------
        fasta :
            an input fasta file from a list of SeqRecord objects of originally string the amino acid sequences

        Raises
        -------
        TypeError
        """
        records = []
        for index, seq in lists.itertuples():
            record = SeqRecord(Seq(seq), id=str(index))
            records.append(record)
        
        with open(f"{inputname}.fasta", "w") as file:
                SeqIO.write(records, file, "fasta")
        return file

    def run_hmmer(DB_PATH: str, input_filename: str, output_filename: str, num_cpu: int = 4) -> None:
        """
        Runs HMMER search against the Pfam database using the hmmscan command
    
        Parameters:
        ------------
        DB_PATH : str
            the path of the HMMER/pfam database
        input_filename : str
            a name for the input fasta file (don't forget extension!)
        output_filename : str
            a name for the domtblout output file (don't forget extension!)
        num_cpu : int

    
        Returns:
        ------------
        domtblout : 
            an output domtblout file of the HMMER/pfam results

        Raises
        -------
        TypeError
        """
        with open(output_filename, "w", encoding="utf-8") as file:
            subprocess.run(["hmmscan", "--cpu", str(num_cpu), "--domtblout", output_filename, DB_PATH, input_filename], stdout=file)
        return file

    # Define a function to find the best hit (hsps: high-scoring segment pair)
    def find_best_hit(query_id, hits):
        """"
        Given a list of HMMER hits, returns the best hit based on e-value.
    
        Parameters:
        ------------
        query_id : str 
            The ID of the query sequence.
        hits : list 
            List of HMMER hits to search through.
    
        Returns:
        ------------
        dict :
            A dictionary containing information about the best hit, including family ID, e-value, 
            sequence length, bit score, high-scoring pair bias, domain e-value, start position, and end position.
            If no hit is found, returns None.
    
        Raises
        -------
    
        """
        best_hit = None
        for hit in hits:
            for hsp in hit.hsps:
                family_id = hit.id
                evalue = hit.evalue
                seqlen = hit.seq_len
                bitscore = hsp.bitscore
                devalue = hsp.evalue
                bias = hsp.bias
                start = hsp.query_start
                end = hsp.query_end
                if best_hit is None or evalue < best_hit["evalue"]:
                    best_hit = {"query_id": query_id, "family_id": family_id, "evalue": evalue,"seqlen": seqlen ,"bias": bias, "devalue": devalue, "bitscore": bitscore, "start": start, "end": end}
        return best_hit

    # Parse the HMMER results file using HmmerIO
    # Store the counts of hits for each family ID
    def parse_hmmer(output_name: str):
        """
        parses HMMER output result

        Parameters:
        ------------
        output_name : str
            a name for the domtblout output file (don't forget extension!)

    
        Returns:
        ------------
        file : 
            an input fasta file from a list of SeqRecord objects of originally string the amino acid sequences

        Raises
        -------
        TypeError
        """
        # Parse the HMMER results file using HmmerIO
        with open(output_name, "r") as results_file:
            qresults = list(SearchIO.parse(results_file, 'hmmscan3-domtab'))

        # Store the counts of hits for each family ID
        hit_counts = defaultdict(int)
        for qresult in qresults:
            query_id = qresult.id
            best_hit = find_best_hit(query_id, qresult.hits)
            if best_hit is not None:
                family_id = best_hit["family_id"]
                hit_counts[family_id] += 1
                print(f"Query sequence ID: {query_id}")
                print(f"Best Pfam family match: {family_id}")
                print(f"Total E-value: {best_hit['evalue']}")
                print(f"Length of hit sequence or HMM: {best_hit['seqlen']}")
                print(f"High-scoring pair bias: {best_hit['bias']}")
                print(f"Domain E-vlaue: {best_hit['devalue']}")
                print(f"Bit score: {best_hit['bitscore']}")
                print(f"Start position: {best_hit['start']}")
                print(f"End position: {best_hit['end']}")
                print("----")
            else:
                print("No Pfam family match found for query sequence ID:", query_id)
    
    ### wrapper execution
    # User should define an amino acid sequence as Pandas series
    
    # generate meso and thermo files
    read_seq(seq, input_filename)

    # place files into HMMER/pfam
    run_hmmer(pfam_path, input_filename_with_ext, output_filename_with_ext, cpu)

    # parse outputs of HMMER/pfam
    parse_hmmer(output_filename_with_ext)



### Data prep and processing

# reading the data
protein_database = pd.read_csv(ID_DB_PATH, index_col=0)

# separating the meso and thermo
meso_seq_db = protein_database.loc[protein_database["identity"] == "m"]
thermo_seq_db = protein_database.loc[protein_database["identity"] == "t"]

## processing meso to be suitable for HMMER
meso_seq_pre = meso_seq_db[["protein_seq", "protein_int_index"]]
meso_seq_list = meso_seq_pre.set_index("protein_int_index").iloc[:15]
meso_seq_list.index.name = None


## processing meso to be suitable for HMMER
thermo_seq_pre = thermo_seq_db[["protein_seq", "protein_int_index"]]
thermo_seq_list = thermo_seq_pre.set_index("protein_int_index").iloc[:15]
thermo_seq_list.index.name = None


#### Execution
hmmer_wrapper(meso_seq_list, "meso_input", PFAM_PATH, "meso_input.fasta", "meso_output.domtblout", 4)
hmmer_wrapper(thermo_seq_list, "thermo_input", PFAM_PATH, "thermo_input.fasta", "thermo_output.domtblout", 4)
