"""
Using 'OMA-pp_test.db' as a test database, this script will attempt to add a few new entries to the database.
The main goal here is to make sure that we expand the negative entries in the database, so that we can have  
a larger negative space to work with.

Personal notes:
----------------
The pairs table in OMA-pp_test.db (~61k entries) has 19 columns with the following schema:
m_protein_seq, t_protein_seq, meso_alphafold_id, thermo_alphafold_id, meso_pid, thermo_pid
bit_score, local_gap_compressed_percent_id, scaled_local_query_percent_id, scaled_local_symmetric_percent_id
query_align_len, query_align_cov, subject_align_len, subject_align_cov, m_ogt, t_ogt, ogt_difference
m_protein_len, t_protein_len

The protein_pairs table in l2t_500k.db (~200k entries) has 17 columns with the following schema:
thermo_pid, meso_pid, local_gap_compressed_percent_id, scaled_local_query_percent_id, scaled_local_symmetric_percent_id
local_E_value, query_align_start, query_align_end, subject_align_start, subject_align_end, query_align_len, query_align_cov
subject_align_len, subject_align_cov, bit_score, thermo_taxid, meso_taxid
"""
# system dependencies
import os

# library dependencies
import duckdb as ddb

# local dependencies
import pairpro.utils as pp_utils

####################
### PATHS & VARS ###
####################

# Paths

# venv variables
if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
else:
    LOGLEVEL = 'DEBUG' # change to INFO for production
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'



if __name__ == "__main__":
    # initialize logger
    logger = pp_utils.start_logger_if_necessary(
        LOGNAME, LOGFILE, LOGLEVEL, filemode='w'
    )
    logger.info(f"Starting {__file__}...")

    # connect to database
    conn_oma_pp = ddb.connect('OMA-pp_test.db', read_only=False)
    # conn_l2t_full = ddb.connect('learn2therm.ddb', read_only=True)
    conn_l2t_subset = ddb.connect('l2t_500k.db', read_only=True)
    logger.info(f"Successfully connected to the databases.")

    # Step 1: create 'analysis' table with the overlapping schema
    create_table_query = """CREATE OR REPLACE TABLE analysis (
        m_protein_seq VARCHAR,
        t_protein_seq VARCHAR,
        thermo_pid VARCHAR,
        meso_pid VARCHAR,
        local_gap_compressed_percent_id DOUBLE,
        scaled_local_query_percent_id DOUBLE,
        scaled_local_symmetric_percent_id DOUBLE,
        query_align_len BIGINT,
        query_align_cov DOUBLE,
        subject_align_len BIGINT,
        subject_align_cov DOUBLE,
        bit_score DOUBLE,
        true_pair BOOLEAN
        );
        """
    conn_oma_pp.execute(create_table_query)
    logger.info("Successfully created the 'analysis' table in the OMA-pp overlap database.")

    # Step 2: Insert "True" pairs into the 'analysis' table
    # The "True" pairs are the pairs that are in the pairs table in OMA-pp_test.db
    insert_true_pairs_query = """
    INSERT INTO analysis
        (m_protein_seq, t_protein_seq, thermo_pid, meso_pid, local_gap_compressed_percent_id,
        scaled_local_query_percent_id, scaled_local_symmetric_percent_id, 
        query_align_len, query_align_cov, subject_align_len, subject_align_cov, bit_score, true_pair)
    SELECT m_protein_seq, t_protein_seq, thermo_pid, meso_pid, local_gap_compressed_percent_id,
    scaled_local_query_percent_id, scaled_local_symmetric_percent_id, 
    query_align_len, query_align_cov, subject_align_len, subject_align_cov, bit_score, TRUE
    FROM pairs;
    """
    conn_oma_pp.execute(insert_true_pairs_query)

    # Commit changes to the database
    conn_oma_pp.commit()
    logger.info("Successfully migrated the 'pairs' table (true pairs) into the 'analysis' table. Committed changes.")

    # Step 3: Insert "False" pairs into the 'analysis' table. Create 50/50 ratio of true/false pairs.
    conn_oma_pp.execute("ATTACH DATABASE 'l2t_500k.db' AS l2t_subset;")
    insert_negative_pairs_query = """
    INSERT INTO analysis (
        m_protein_seq,
        t_protein_seq,
        thermo_pid,
        meso_pid,
        local_gap_compressed_percent_id,
        scaled_local_query_percent_id,
        scaled_local_symmetric_percent_id,
        query_align_len,
        query_align_cov,
        subject_align_len,
        subject_align_cov,
        bit_score,
        true_pair
    )
    SELECT
        proteins_m.protein_seq,
        proteins_t.protein_seq,
        l2t.thermo_pid,
        l2t.meso_pid,
        l2t.local_gap_compressed_percent_id,
        l2t.scaled_local_query_percent_id,
        l2t.scaled_local_symmetric_percent_id,
        l2t.query_align_len,
        l2t.query_align_cov,
        l2t.subject_align_len,
        l2t.subject_align_cov,
        l2t.bit_score,
        FALSE
    FROM l2t_subset.protein_pairs AS l2t
    JOIN l2t_subset.proteins AS proteins_m ON l2t.thermo_pid = proteins_m.pid
    JOIN l2t_subset.proteins AS proteins_t ON l2t.meso_pid = proteins_t.pid
    WHERE NOT EXISTS (
        SELECT 1 
        FROM analysis 
        WHERE analysis.meso_pid = l2t.meso_pid 
        AND analysis.thermo_pid = l2t.thermo_pid
    )
    LIMIT 63306;
    """
    conn_oma_pp.execute(insert_negative_pairs_query)
    conn_oma_pp.commit()

    # Step 4: data verification
    # check the total number of entries in the 'analysis' table
    row_count = conn_oma_pp.execute("SELECT COUNT(*) FROM analysis;").fetchone()[0]
    logger.info(f"Total number of entries in the 'analysis' table: {row_count}")

    # check the number of true/false pairs
    true_pair_count = conn_oma_pp.execute("SELECT COUNT(*) FROM analysis WHERE true_pair = TRUE;").fetchone()[0]
    negative_pair_count = conn_oma_pp.execute("SELECT COUNT(*) FROM analysis WHERE true_pair = FALSE;").fetchone()[0]
    logger.info(f"Number of 'true' (true_pairs = TRUE) pairs: {true_pair_count}")
    logger.info(f"Number of 'false' (false_pairs = FALSE) pairs: {negative_pair_count}")

    sample_rows = conn_oma_pp.execute("SELECT * FROM analysis LIMIT 5;").fetchall()
    for row in sample_rows:
        logger.info(row)

    logger.info("Successfully created the 'analysis' table with true/false pairs. Committed changes.")
    # check for duplicates in the 'analysis' table
    duplicates_pairs_query = """
    SELECT meso_pid, thermo_pid, COUNT(*) as count_star
    FROM analysis
    GROUP BY meso_pid, thermo_pid
    HAVING COUNT(*) > 1;
    """

    duplicate_pairs = conn_oma_pp.execute(duplicates_pairs_query).fetchall()

    if len(duplicate_pairs) > 0:
        logger.warning(f"Found {len(duplicate_pairs)} duplicate pairs in the 'analysis' table.")
        for pair in duplicate_pairs:
            logger.warning(f"Duplicate pair: meso_pid={pair[0]}, thermo_pid={pair[1]}, count={pair[2]}")
    else:
        logger.info("No duplicate pairs found in the 'analysis' table!")


    conn_oma_pp.close()
    conn_l2t_subset.close()
    logger.info("Successfully closed the database connections.")