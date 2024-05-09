"""
TODO
"""
# system dependencies
import os
import sys
import logging

# library dependencies
import duckdb as ddb
import pandas as pd
from timeit import default_timer as timer

# local dependencies
import pairpro.utils as pp_utils
import pairpro.user_blast as pp_up

####################
### PATHS & VARS ###
####################
# paths
OMA_DIR = "./tmp/oma.db"

# venv variables
if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
else:
    LOGLEVEL = 'DEBUG' # change to INFO for production
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'

# CLI arguments for sbatch
num_cpus = int(sys.argv[1])

if __name__ == "__main__":
    # initialize logger
    logger = pp_utils.start_logger_if_necessary(
        LOGNAME, LOGFILE, LOGLEVEL, filemode='w')
    logger.info("Starting script. Logging to %s", LOGFILE)

    # connect to database
    con = ddb.connect(OMA_DIR, read_only=False)
    logger.info("Connected to database %s", OMA_DIR)

    start_time = timer()
    logger.info(f'Starting sample to generate nuprok pairs; 12M pairs')
    pp_utils.nuprok_sample(con=con, size=12000000)
    logger.info(f'Completed nuprok pairs: {timer()-start_time} seconds')

    logger.info(f"Proceeding to process and clean the table")

    # create 'clean' nuprok pairs
    start_time = timer()
    con.execute("""CREATE OR REPLACE TABLE clean_nuprok_pairs 
                AS SELECT * FROM nuprok_pairs 
                WHERE protein1_uniprot_id <> protein2_uniprot_id AND protein1_sequence <> protein2_sequence 
                AND protein1_uniprot_id IS NOT NULL 
                AND protein2_uniprot_id IS NOT NULL 
                AND protein1_sequence IS NOT NULL 
                AND protein2_sequence IS NOT NULL""")
    logger.info(f'Completed processing nuprok pairs: {timer()-start_time} seconds')
    logger.info(f'clean_nuprok_pairs table created. {con.execute("SELECT COUNT(*) FROM clean_nuprok_pairs").fetchdf()} rows, and has the following columns: {con.execute("DESCRIBE clean_nuprok_pairs").fetchall()}')

    # create 'bad' sample from clean_nuprok_pairs
    start_time = timer()
    combined_pairs_query = """
    WITH shuffled_protein1 AS (
        SELECT
            row_number() OVER (ORDER BY random()) AS rn,
            protein1_uniprot_id, protein1_sequence
        FROM clean_nuprok_pairs
        USING SAMPLE 100%
    ),
    shuffled_protein2 AS (
        SELECT
            row_number() OVER (ORDER BY random()) AS rn,
            protein2_uniprot_id, protein2_sequence
        FROM clean_nuprok_pairs
        USING SAMPLE 100%
    ),
    bad_pairs AS (
        SELECT
            'bad_' || sp1.rn AS pair_id,
            sp1.protein1_uniprot_id, sp2.protein2_uniprot_id,
            sp1.protein1_sequence, sp2.protein2_sequence
        FROM shuffled_protein1 sp1
        JOIN shuffled_protein2 sp2 ON sp1.rn = sp2.rn
    ),
    clean_pairs AS (
        SELECT
            'clean_' || row_number() OVER () AS pair_id,
            protein1_uniprot_id, protein2_uniprot_id,
            protein1_sequence, protein2_sequence
        FROM clean_nuprok_pairs
    )
    SELECT * FROM clean_pairs WHERE LENGTH(protein1_sequence) < 250 AND LENGTH(protein2_sequence) < 250
    UNION ALL
    SELECT * FROM bad_pairs WHERE LENGTH(protein1_sequence) < 250 AND LENGTH(protein2_sequence) < 250
    """

    # execute the query to combine clean and bad pairs
    con.execute(f"""CREATE OR REPLACE TABLE combined_pairs AS ({combined_pairs_query})""")
    logger.info(f'Completed creating combined_pairs: {timer()-start_time} seconds')


    logger.info(f"Proceeding to align nuprok pairs. Making into dataframe to BLAST.")
    df = con.execute("SELECT * FROM combined_pairs").fetchdf()

    # align combined pairs
    start_time = timer()
    module_path=os.getcwd()
    logger.debug(f"CHECK-> module path: {module_path}")

    logger.info("Starting to align nuprok pairs")
    df_align = pp_up.make_blast_df(df, cpus=num_cpus, path='./data/OMA/samples/oma_sample_12m.db')
    time_delta = timer() - start_time
    logger.info(f'Alignment time: {time_delta} seconds using {num_cpus} cpus')
    
    
    logger.info("Completed database operations!!")

    # save aligned dataframe
    logger.info("Saving aligned dataframe...")
    df_align[['protein1_alphafold_id', 'protein2_alphafold_id']] = df[['protein1_uniprot_id', 'protein2_uniprot_id']]
    df_align.to_csv('./data/OMA/samples/oma_sample_12m_aligned.csv', index=False)
    logger.info("Completed alignment. Saved to ./data/OMA/samples/oma_sample_12m_aligned.csv")