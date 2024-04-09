"""
An analysis script
"""
# system dependencies
import os

# library dependencies
import duckdb as ddb
import pandas as pd

# local dependencies
import pairpro.utils as pp_utils

####################
### PATHS & VARS ###
####################


# venv variables
if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
else:
    LOGLEVEL = 'DEBUG' # change to INFO for production
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'

####################
### MAIN PROGRAM ###
####################


if __name__ == "__main__":
    # initialize logger
    logger = pp_utils.start_logger_if_necessary(
        LOGNAME, LOGFILE, LOGLEVEL, filemode='w'
    )
    logger.info(f"Running {__file__}...")

    # Step 0: connect to the databases
    conn_pairpro_50k = ddb.connect('./tmp/pairpro.db', read_only=True)
    # conn_pairpro_500k = ddb.connect('./tmp/pairpro_2.db', read_only=True) ## This is 226k pairs db
    conn_oma_pp = ddb.connect('OMA-pp_test.db', read_only=True)
    logger.info("Successfully connected to the databases.")

    # 1. Count the total number of protein pairs in each of the databases
    pairpro_unique_protein_pairs = conn_pairpro_50k.execute("SELECT COUNT(*) FROM pairpro.pairpro.final").fetchone()[0]
    # pairpro_unique_protein_pairs_500k = conn_pairpro_500k.execute("SELECT COUNT(*) FROM pairpro.pairpro.final").fetchone()[0]
    oma_pairs_count = conn_oma_pp.execute("SELECT COUNT(*) FROM analysis;").fetchone()[0]

    # 2. Count the number of unique proteins in each of the databases
    pairpro_unique_protein_count = conn_pairpro_50k.execute("SELECT COUNT(DISTINCT pid) FROM pairpro.pairpro.proteins;").fetchone()[0]
    # pairpro_unique_protein_count_500k = conn_pairpro_500k.execute("SELECT COUNT(DISTINCT pid) FROM pairpro.pairpro.proteins;").fetchone()[0]
    oma_unique_protein_count = conn_oma_pp.execute("SELECT COUNT(DISTINCT meso_pid) + COUNT(DISTINCT thermo_pid) FROM analysis;").fetchone()[0]


    logger.info(f"Pairpro_50k - total number of unique protein pairs: {pairpro_unique_protein_pairs}")
    logger.info(f"OMA-pp_test - total number of unique protein pairs: {oma_pairs_count}")
    logger.info("----")
    logger.info(f"Pairpro_50k - total number of unique proteins: {pairpro_unique_protein_count}")
    logger.info(f"OMA-pp_test - total number of unique proteins: {oma_unique_protein_count}")


    # Count the number of true pairs and false pairs in OMA-pp_test.db
    true_pairs_count = conn_oma_pp.execute("SELECT COUNT(*) FROM analysis WHERE true_pair = 1;").fetchone()[0]
    false_pairs_count = conn_oma_pp.execute("SELECT COUNT(*) FROM analysis WHERE true_pair = 0;").fetchone()[0]

    logger.info(f"OMA-pp_test - True pairs: {true_pairs_count}")
    logger.info(f"OMA-pp_test - False pairs: {false_pairs_count}")

    ####### Overlap analysis #######
    logger.info("----")
    logger.info("Overlap analysis...")

    # I got a duckdb binder error, so I will attach the database
    conn_oma_pp.execute("ATTACH DATABASE './tmp/pairpro.db' AS pairpro_50k;")

    # a. Overlapping Protein Pairs
    overlap_pairs_query = """
    SELECT COUNT(*)
    FROM (SELECT meso_pid, thermo_pid FROM pairpro_50k.pairpro.pairpro.final)
    INTERSECT
    (SELECT meso_pid, thermo_pid FROM analysis WHERE true_pair = 1);
    """

    overlap_pairs_count = conn_oma_pp.execute(overlap_pairs_query).fetchone()[0]

    # b. Overlapping Unique Proteins
    overlap_proteins_query = """
    SELECT COUNT(*)
    FROM (SELECT DISTINCT pid FROM pairpro_50k.pairpro.pairpro.proteins)
    INTERSECT
    (SELECT DISTINCT meso_pid FROM analysis UNION SELECT DISTINCT thermo_pid FROM analysis);
    """

    overlap_proteins_count = conn_oma_pp.execute(overlap_proteins_query).fetchone()[0]

    logger.info(f"Number of overlapping protein pairs between Pairpro_50k and OMA-pp_test (true pair): {overlap_pairs_count}")
    logger.info(f"Number of overlapping unique proteins between Pairpro_50k and OMA-pp_test: {overlap_proteins_count}")
