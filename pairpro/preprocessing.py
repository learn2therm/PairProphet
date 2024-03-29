'''
This package builds the PairProphet database from learn2thermDB.

Functions:
    connect_df: Establishes connection to DuckDB database using local or
                remote input path. Reports time to connection.

    build_pairpro: Constructs pairprophet database from learn2therm database.
'''

import time
import duckdb
import os


def connect_db(path: str, empty=False):
    '''
    Runs duckdb.connect() function on database path. Returns a
    duckdb.DuckDBPyConnection object and prints execution time.

    Args:
        path (str): Path to DuckDB database file containing learn2therm.

    Returns:
        con (duckdb.DuckDBPyConnection): A DuckDB connection object linking
                                         script to learn2therm database.

    Raises:
        AttributeError: Input database contains no tables.
    '''
    s_time = time.time()

    print('Connecting to database...')
    con = duckdb.connect(path)

    if empty is False:
        tables = con.execute("""SELECT TABLE_NAME
                                FROM INFORMATION_SCHEMA.TABLES
                                WHERE TABLE_TYPE='BASE TABLE'""").df()
        if tables.shape[0] < 1:
            raise AttributeError('Input database is empty.')
    else:
        tables = []

    e_time = time.time()
    elapsed_time = e_time - s_time
    print(f'Connection established! Execution time: {elapsed_time} seconds')
    return con, tables


def build_pairpro(con, out_db_path, min_ogt_diff: int = 20,
                  min_16s: int = 1300):
    '''
    Converts learn2therm DuckDB database into a DuckDB database for PairProphet
    by adding filtered and constructed tables. Ensure at lease 20 GB of free
    disk space and 30 GB of system memory are available before running on the
    full database.

    Args:
        con (duckdb.DuckDBPyConnection): DuckDB connection object. Links script
                                         to DuckDB SQL database.
        out_db_path (str): Path to PairProphet output database file.
        min_ogt_diff (int): Cutoff for minimum difference in optimal growth
                            temperature between thermophile and mesophile
                            pairs. Default 20 deg C.
        min_16s (int): Cutoff for minimum 16S read length for taxa. Default
                       1300 bp. Filters out organisms with poor or incomplete
                       16S sequencing.

    Returns:
        None. Database object is modified in place.

    Raises:
        ValueError: Optimal growth temperature difference must be positive.
        ValueError: Minimum 16S sequence read is 1 bp.
        AttributeError: Database must be in the learn2therm format.
    '''
    if min_ogt_diff < 0:
        raise ValueError("""Optimal growth temperature difference must be
                         positive.""")

    if min_16s < 1:
        raise ValueError('16S must have at least 1 bp read.')

    tables = con.execute("""SELECT TABLE_NAME
                            FROM INFORMATION_SCHEMA.TABLES
                            WHERE TABLE_TYPE='BASE TABLE'""").df()

    # Check if proper tables exist in database. If they do not, raise an error.
    if (item in tables for item in ['proteins', 'protein_pairs', 'taxa',
                                    'taxa_pairs']):
        pass

    else:
        raise AttributeError('Database is not formatted for learn2therm.')

    s_time = time.time()
    print('Constructing pairpro_taxa_pairs...')

    # Builds PairProphet taxa pair table using paired taxa from learn2therm
    taxa_pairs_cmd = """CREATE OR REPLACE TEMP TABLE pairpro_taxa_pairs AS
                        SELECT *
                        FROM taxa_pairs
                        INNER JOIN taxa_pairs_lab
                        ON (taxa_pairs.__index_level_0__ =
                        taxa_pairs_lab.__index_level_0__)
                        WHERE taxa_pairs_lab.is_pair = True"""
    con.execute(taxa_pairs_cmd)

    e_time = time.time()
    elapsed_time = e_time - s_time
    print(f"""Finished constructing pairpro_taxa_pairs. Execution time:
          {elapsed_time} seconds""")
    print('Constructing pairpro_taxa...')

    # Builds PairProphet taxa table using only paired taxa from learn2therm.
    taxa_cmd = """CREATE OR REPLACE TEMP TABLE pairpro_taxa AS
                  SELECT *
                  FROM taxa
                  WHERE taxid IN
                  (SELECT DISTINCT query_id FROM pairpro_taxa_pairs)
                  OR taxid IN
                  (SELECT DISTINCT subject_id FROM pairpro_taxa_pairs)"""
    con.execute(taxa_cmd)

    e_time2 = time.time()
    elapsed_time = e_time2 - e_time
    print(f"""Finished constructing pairpro_taxa. Execution time:
          {elapsed_time} seconds""")
    print('Filtering on ogt and 16S sequence parameters...')

    # Builds PairProphet table containing taxa pairs and their associated
    # optimal growth temperatures (ogt). Excludes 16S sequences and ogt
    # difference below cutoff values from function input.
    ogt_pairs_cmd = f"""CREATE OR REPLACE TEMP TABLE pairpro_ogt_taxa_pairs AS
                        SELECT pairpro_taxa_pairs.*,
                        taxa_m.temperature AS meso_ogt,
                        taxa_t.temperature AS thermo_ogt,
                        taxa_t.temperature - taxa_m.temperature AS ogt_diff,
                        taxa_m."16s_len" AS meso_16s_len,
                        taxa_t."16s_len" AS thermo_16s_len
                        FROM pairpro_taxa_pairs
                        JOIN pairpro_taxa AS taxa_m
                        ON (pairpro_taxa_pairs.subject_id = taxa_m.taxid)
                        JOIN pairpro_taxa AS taxa_t
                        ON (pairpro_taxa_pairs.query_id = taxa_t.taxid)
                        WHERE ogt_diff >= {min_ogt_diff}
                        AND meso_16s_len >= {min_16s}
                        AND thermo_16s_len >= {min_16s}"""
    con.execute(ogt_pairs_cmd)

    e_time3 = time.time()
    elapsed_time = e_time3 - e_time2
    print(f'Finished filtering. Execution time: {elapsed_time} seconds')
    print('Constructing pairpro_protein_pairs...')

    # Builds PairProphet table containing protein pairs
    protein_pair_cmd = """CREATE OR REPLACE TEMP TABLE pairpro_protein_pairs AS
                          SELECT protein_pairs.*,
                          otp.meso_ogt AS m_ogt,
                          otp.thermo_ogt AS t_ogt,
                          otp.ogt_diff AS ogt_difference
                          FROM protein_pairs
                          INNER JOIN pairpro_ogt_taxa_pairs AS otp
                          ON (protein_pairs.thermo_taxid = otp.query_id)
                          AND (protein_pairs.meso_taxid = otp.subject_id)"""
    con.execute(protein_pair_cmd)

    e_time4 = time.time()
    elapsed_time = e_time4 - e_time3
    print(f"""Finished constructing pairpro_protein_pairs. Execution time:
          {elapsed_time} seconds""")
    print('Constructing pairpro_proteins...')

    # Builds PairProphet table containing proteins that belong to taxa from
    # pairpro_taxa_pairs.
    prot_filt_cmd = """CREATE OR REPLACE TEMP TABLE pairpro_proteins AS SELECT *
                       FROM proteins
                       WHERE pid IN (SELECT DISTINCT meso_pid
                                     FROM pairpro_protein_pairs)
                       OR
                       pid IN (SELECT DISTINCT thermo_pid
                               FROM pairpro_protein_pairs)
                    """
    con.execute(prot_filt_cmd)

    e_time5 = time.time()
    elapsed_time = e_time5 - e_time4
    print(f"""Finished constructing pairpro_proteins. Execution time:
          {elapsed_time} seconds""")
    print('Constructing final dataset...')

    # Builds final PairProphet data table for downstream sampling.
    big_table_cmd = """CREATE OR REPLACE TEMP TABLE pairpro_final AS
                       SELECT pairpro_protein_pairs.*,
                       proteins_m.protein_seq AS m_protein_seq,
                       proteins_t.protein_seq AS t_protein_seq,
                       proteins_m.pdb_id AS meso_pdb,
                       proteins_t.pdb_id AS thermo_pdb
                       FROM pairpro_protein_pairs
                       JOIN pairpro_proteins AS proteins_m
                       ON (pairpro_protein_pairs.meso_pid = proteins_m.pid)
                       JOIN pairpro_proteins AS proteins_t
                       ON (pairpro_protein_pairs.thermo_pid =
                           proteins_t.pid)"""
    con.execute(big_table_cmd)

    if os.path.exists(out_db_path):
        filename, ext = os.path.splitext(out_db_path)
        counter = 1
        while os.path.exists(f'{filename}_{counter}{ext}'):
            counter += 1
        out_db_path = f'{filename}_{counter}{ext}'
        filename = f'{filename.split("/")[2]}_{counter}'
    else:
        filename = os.path.splitext(out_db_path)[0].split("/")[2]

    print(f'Transferring data to new database {out_db_path}')
    con.execute(f"""ATTACH '{out_db_path}' AS out_db""")
    con.execute("""CREATE SCHEMA out_db.pairpro""")
    con.execute("""CREATE OR REPLACE TABLE out_db.pairpro.final AS
                   SELECT * FROM main.pairpro_final""")
    con.execute("""CREATE OR REPLACE TABLE out_db.pairpro.proteins AS
                   SELECT * FROM main.pairpro_proteins""")
    con.execute("""DETACH out_db""")

    print('Finishing up...')
    con.commit()
    con.close()

    con2, _ = connect_db(out_db_path)

    # Add pair IDs to final table
    con2.execute(f"""CREATE TEMP TABLE pair_ids AS
                     SELECT ROW_NUMBER() OVER(ORDER BY meso_pid, thermo_pid)
                     AS pair_id, meso_pid, thermo_pid
                     FROM {filename}.pairpro.final""")
    con2.execute(f"""ALTER TABLE {filename}.pairpro.final ADD pair_id int""")
    con2.execute(f"""UPDATE {filename}.pairpro.final AS f
    SET pair_id = pair_ids.pair_id::int
    FROM pair_ids
    WHERE pair_ids.meso_pid = f.meso_pid
    AND pair_ids.thermo_pid = f.thermo_pid
    """)

    et_final = time.time()
    elapsed_time = et_final - e_time5
    print(f'Finished. Total execution time: {elapsed_time} seconds')

    return con2, filename
