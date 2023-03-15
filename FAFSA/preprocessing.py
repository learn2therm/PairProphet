'''
This package builds the FAFSA database from learn2therm.

Functions:
    connect_df: Establishes connection to DuckDB database using local or
                remote input path. Reports time to connection.

    build_fafsa: Constructs fafsa database from learn2therm database.

    sankey_plots: Optional function run by calling build_fafsa with plots = True.
                  Generates Sankey plots showing the fate of samples as they pass from
                  learn2therm tables to FAFSA tables.
'''

import os
import time

import numpy as np

import duckdb

# Dependencies for Sankey plots.
import plotly.graph_objects as go
import plotly.io as pio
pio.kaleido.scope.default_format = "png"


def connect_db(path: str):
    '''
    Runs duckdb.connect() function on database path. Returns a
    duckdb.DuckDBPyConnection object and prints execution time.

    Args:
        path (str): Path to DuckDB database file containing learn2therm.

    Returns:
        con (duckdb.DuckDBPyConnection): A DuckDB connection object linking script to
        learn2therm database.

    Raises:
        VersionError: DuckDB installation is not one of 0.6.0 or 0.6.1.
        AttributeError: Input database contains no tables.
    '''
    s_time = time.time()

    version = duckdb.__version__
    print(f'Version {version} detected.')

    # Checks for compatible installation of duckdb.
    if not version in ['0.6.0', '0.6.1']:
        raise VersionError("""learn2therm was generated using DuckDB storage version 39. It is only
                           compatible with duckdb versions 0.6.0 and 0.6.1. Please check your
                           installation. Refer to https://duckdb.org/internals/storage.html for more
                           details.""")

        e_time = time.time()
        elapsed_time = e_time - s_time
        print(f'Finished with VersionError. Execution time: {elapsed_time} seconds')

    print('Connecting to database...')
    con = duckdb.connect(path)

    tables = con.execute("""SELECT TABLE_NAME
                            FROM INFORMATION_SCHEMA.TABLES
                            WHERE TABLE_TYPE='BASE TABLE'""").df()

    if tables.shape[0] < 1:
        raise AttributeError('Input database is empty.')

    e_time = time.time()
    elapsed_time = e_time - s_time
    print(f'Connection established! Execution time: {elapsed_time} seconds')
    return con


def build_fafsa(con, min_ogt_diff: int = 20, min_16s: int = 1300,
                    plots: bool = False):
    '''
    Converts learn2therm DuckDB database into a DuckDB database for FAFSA by adding filtered and
    constructed tables. Ensure at lease 100 GB of free disk space and 30 GB of system memory are 
    available before running on the full database.

    Args:
        con (duckdb.DuckDBPyConnection): DuckDB connection object. Links script to DuckDB SQL 
        database.
        min_ogt_diff (int): Cutoff for minimum difference in optimal growth temperature between 
        thermophile and mesophile pairs. Default 20 deg C.
        min_16s (int): Cutoff for minimum 16S read length for taxa. Default 1300 bp. Filters out 
        organisms with poor or incomplete 16S sequencing.
        plots (bool): Boolean to determine whether the user wants Sankey plots diagramming the fate
        of learn2therm samples during FAFSA construction to be saved in ./plots.

    Returns:
        None. Database object is modified in place.

    Raises:
        ValueError: Optimal growth temperature difference must be positive.
        ValueError: Minimum 16S sequence read is 1 bp.
        AttributeError: Database must be in the learn2therm format.
    '''

    if min_ogt_diff < 0:
        raise ValueError('Optimal growth temperature difference must be positive.')

    if min_16s < 1:
        raise ValueError('16S must have at least 1 bp read.')

    tables = con.execute("""SELECT TABLE_NAME
                            FROM INFORMATION_SCHEMA.TABLES
                            WHERE TABLE_TYPE='BASE TABLE'""").df()

    # Check if proper tables exist in database. If they do not, raise an error.
    if (item in tables for item in ['proteins', 'protein_pairs', 'taxa', 'taxa_pairs']):
        pass

    else:
        raise AttributeError('Database is not formatted for learn2therm.')

    s_time = time.time()
    print('Constructing fafsa_taxa_pairs...')

    # Builds FAFSA taxa pair table using only paired taxa from learn2therm
    taxa_pairs_cmd = """CREATE OR REPLACE TABLE fafsa_taxa_pairs AS
                        SELECT *
                        FROM taxa_pairs
                        WHERE is_pair = True"""
    con.execute(taxa_pairs_cmd)

    e_time = time.time()
    elapsed_time = e_time - s_time
    print(f'Finished constructing fafsa_taxa_pairs. Execution time: {elapsed_time} seconds')
    print('Constructing fafsa_taxa...')

    # Commands to identify all taxa that are implicated in learn2therm pairs.
    meso_cmd = """SELECT DISTINCT meso_index
                  FROM taxa_pairs
                  WHERE is_pair = True"""
    thermo_cmd = """SELECT DISTINCT thermo_index
                    FROM taxa_pairs
                    WHERE is_pair = True"""

    useful_thermo = con.execute(thermo_cmd).df()
    useful_meso = con.execute(meso_cmd).df()

    # Generates tuple object containing all relevant taxa
    useful_taxa = tuple(list(useful_meso['meso_index']) + list(useful_thermo['thermo_index']))

    # Builds FAFSA taxa table using only paired taxa from learn2therm.
    taxa_cmd = f"""CREATE OR REPLACE TABLE fafsa_taxa AS
                   SELECT *
                   FROM taxa
                   WHERE taxa_index IN {useful_taxa}"""
    con.execute(taxa_cmd)

    e_time2 = time.time()
    elapsed_time = e_time2 - e_time
    print(f'Finished constructing fafsa_taxa. Execution time: {elapsed_time} seconds')
    print('Filtering on ogt and 16S sequence parameters...')

    # Builds FAFSA table containing taxa pairs and their associated optimal growth temperatures
    # (ogt). Excludes 16S sequences and ogt difference below cutoff values from function input.
    ogt_pairs_cmd = f"""CREATE OR REPLACE TABLE fafsa_ogt_taxa_pairs AS SELECT fafsa_taxa_pairs.*,
                        taxa_m.ogt AS meso_ogt,
                        taxa_t.ogt AS thermo_ogt,
                        taxa_t.ogt - taxa_m.ogt AS ogt_diff,
                        taxa_m.len_16s AS meso_16s_len,
                        taxa_t.len_16s AS thermo_16s_len
                        FROM fafsa_taxa_pairs
                        JOIN fafsa_taxa AS taxa_m ON (fafsa_taxa_pairs.meso_index = taxa_m.taxa_index)
                        JOIN fafsa_taxa AS taxa_t ON (fafsa_taxa_pairs.thermo_index = taxa_t.taxa_index)
                        WHERE ogt_diff >= {min_ogt_diff}
                        AND meso_16s_len >= {min_16s}
                        AND thermo_16s_len >= {min_16s}"""
    con.execute(ogt_pairs_cmd)

    e_time3 = time.time()
    elapsed_time = e_time3 - e_time2
    print(f'Finished filtering. Execution time: {elapsed_time} seconds')
    print('Constructing fafsa_protein_pairs...')

    # Builds FAFSA table containing protein pairs
    protein_pair_cmd = """CREATE OR REPLACE TABLE fafsa_protein_pairs AS
                          SELECT protein_pairs.*,
                          otp.local_gap_compressed_percent_id AS local_gap_compressed_percent_id_16s,
                          otp.scaled_local_query_percent_id AS scaled_local_query_percent_id_16s,
                          otp.scaled_local_symmetric_percent_id AS scaled_local_symmetric_percent_id_16s,
                          otp.query_align_cov AS query_align_cov_16s,
                          otp.subject_align_cov AS subject_align_cov_16s,
                          otp.bit_score AS bit_score_16s,
                          otp.meso_ogt AS m_ogt,
                          otp.thermo_ogt AS t_ogt,
                          otp.ogt_diff AS ogt_difference
                          FROM protein_pairs
                          INNER JOIN fafsa_ogt_taxa_pairs AS otp
                          ON (protein_pairs.taxa_pair_index = otp.taxa_pair_index)"""
    con.execute(protein_pair_cmd)

    e_time4 = time.time()
    elapsed_time = e_time4 - e_time3
    print(f'Finished constructing fafsa_protein_pairs. Execution time: {elapsed_time} seconds')
    print('Constructing fafsa_proteins...')

    # Builds FAFSA table containing proteins that belong to taxa from fafsa_taxa_pairs.
    prot_filt_cmd = """CREATE OR REPLACE TABLE fafsa_proteins AS SELECT *
                       FROM proteins
                       WHERE protein_int_index IN (SELECT DISTINCT meso_protein_int_index FROM protein_pairs) OR
                       protein_int_index IN (SELECT DISTINCT thermo_protein_int_index FROM protein_pairs)
                    """
    con.execute(prot_filt_cmd)

    e_time5 = time.time()
    elapsed_time = e_time5 - e_time4
    print(f'Finished constructing fafsa_proteins. Execution time: {elapsed_time} seconds')
    print('Constructing final dataset...')

    # Builds final FAFSA data table for downstream sampling.
    big_table_cmd = """CREATE OR REPLACE TABLE fafsa_final AS
                       SELECT fafsa_protein_pairs.*,
                       proteins_m.protein_seq AS m_protein_seq,
                       proteins_t.protein_seq AS t_protein_seq,
                       proteins_m.protein_desc AS m_protein_desc,
                       proteins_t.protein_desc AS t_protein_desc,
                       proteins_m.protein_len AS m_protein_len,
                       proteins_t.protein_len AS t_protein_len
                       FROM fafsa_protein_pairs
                       JOIN fafsa_proteins AS proteins_m
                       ON (fafsa_protein_pairs.meso_protein_int_index = proteins_m.protein_int_index)
                       JOIN fafsa_proteins AS proteins_t
                       ON (fafsa_protein_pairs.thermo_protein_int_index =
                           proteins_t.protein_int_index)"""
    con.execute(big_table_cmd)

    if plots is True:

        sankey_plots(con, min_ogt_diff)

    else:
        pass

    print('Finishing up...')
    con.commit()
    con.close()

    et_final = time.time()
    elapsed_time = et_final - e_time5
    print(f'Finished. Total execution time: {elapsed_time} seconds')


def sankey_plots(con, min_ogt_diff):
    '''
    Constructs Sankey plots for learn2therm to FAFSA data flow. Saves plots to new or existing
    ./plots folder. Does not show 16S filtering yet. As long as 16S cutoff is default 1300 bp, this
    will not affect the accuracy of Sankey plots since the original learn2therm database was
    filtered on this value during construction.

    Args:
        con (duckdb.DuckDBPyConnection): DuckDB connection object. Links script to DuckDB SQL
                                         database.
        min_ogt_diff (int): Cutoff for minimum difference in optimal growth temperature between
                            thermophile and mesophile pairs. Default 20 deg C.

    Returns:
        None. Plots are automatically saved as png in ./plots

    Raises:
    '''

    if min_ogt_diff < 0:
        raise ValueError('Optimal growth temperature difference must be positive.')

    print('Constructing plots.')

    # Checks if plots directory exists and makes one if it does not.
    newdir = 'plots'
    exist = os.path.exists(os.path.join('.', newdir))

    if not exist:
        os.makedirs(newdir)

    # Color palette for Sankey plots.
    sank_blue = 'rgb(32, 159, 223)'
    sank_red = 'rgb(204, 102, 119)'
    sank_grey = 'rgb(221, 221, 221)'
    sank_purple = 'rgb(118, 111, 159)'
    sank_blue_t = 'rgba(32, 159, 223, 0.5)'
    sank_red_t = 'rgba(204, 102, 119, 0.5)'
    sank_grey_t = 'rgba(221, 221, 221, 0.5)'
    sank_purple_t = 'rgba(118, 111, 159, 0.5)'

    # Parameters for taxa pair Sankey
    size_tp_l2t = int(con.execute("""SELECT COUNT(taxa_pair_index)
                                     FROM taxa_pairs""").df().values)
    size_tp_no_pair = int(con.execute("""SELECT COUNT(taxa_pair_index)
                                         FROM taxa_pairs
                                         WHERE is_pair = False""").df().values)
    size_tp_16s_pair = int(con.execute("""SELECT COUNT(taxa_pair_index)
                                          FROM taxa_pairs
                                          WHERE is_pair = True""").df().values)
    size_tp_fafsa = int(con.execute("""SELECT COUNT(taxa_pair_index)
                                           FROM fafsa_taxa_pairs""").df().values)
    size_tp_small_diff = size_tp_16s_pair - size_tp_fafsa

    perc_tp_no_pair = np.round(100*size_tp_no_pair/size_tp_l2t, 1)
    perc_tp_fafsa = np.round(100*size_tp_fafsa/size_tp_l2t, 1)
    perc_tp_small_diff = np.round(100*size_tp_small_diff/size_tp_l2t, 1)

    # Builds and saves taxa pair Sankey in data folder
    fig1 = go.Figure(data=[go.Sankey(
    arrangement = 'snap',
    node = dict(
      pad = 25,
      thickness = 20,
      line = dict(color = 'black', width = 0.5),
      label = [f'learn2therm {str(np.round(size_tp_l2t/1000, 1))+"k"}', '16S pair',
               f'No pair ({perc_tp_no_pair}%)',
               f"""FAFSA ({str(np.round(size_tp_fafsa/1000, 1))+"k"},
               {perc_tp_fafsa}%)""", f"""< {min_ogt_diff} \N{DEGREE SIGN}C
               diff ({perc_tp_small_diff}%)"""],
      color = [sank_purple, sank_purple, sank_grey, sank_purple, sank_grey]
    ),
    link = dict(
      source = [1, 0, 0, 1],
      target = [4, 1, 2, 3],
      value = [size_tp_small_diff, size_tp_16s_pair, size_tp_no_pair, size_tp_fafsa],
      color = [sank_grey_t, sank_purple_t, sank_grey_t, sank_purple_t]
      ))])

    fig1.update_layout(title_text='Taxa Pairs', font_family = 'Arial', font_size=16)
    fig1.write_image(os.path.join('../data', newdir, 'taxa_pair_sankey.png'), engine = 'kaleido',
                     scale = 6, width = 1280, height = 640)

    # Parameters for taxa Sankey
    size_t_l2t = int(con.execute("""SELECT COUNT(DISTINCT taxa_index)
                                    FROM taxa""").df().values)
    size_t_meso = int(con.execute("""SELECT COUNT(DISTINCT meso_index)
                                     FROM taxa_pairs""").df().values)
    size_t_thermo = int(con.execute("""SELECT COUNT(DISTINCT thermo_index)
                                       FROM taxa_pairs""").df().values)
    size_tm_no_pair = int(con.execute("""SELECT COUNT(taxa_index)
                                         FROM taxa
                                         WHERE taxa_index IN
                                         (SELECT meso_index
                                         FROM taxa_pairs
                                         WHERE is_pair = False)""").df().values)
    size_tt_no_pair = int(con.execute("""SELECT COUNT(taxa_index)
                                         FROM taxa
                                         WHERE taxa_index IN
                                         (SELECT thermo_index
                                         FROM taxa_pairs
                                         WHERE is_pair = False)""").df().values)
    size_tm_16s_pair = int(con.execute("""SELECT COUNT(taxa_index)
                                          FROM taxa
                                          WHERE taxa_index IN
                                          (SELECT meso_index
                                          FROM taxa_pairs
                                          WHERE is_pair = True)""").df().values)
    size_tt_16s_pair = int(con.execute("""SELECT COUNT(taxa_index)
                                          FROM taxa
                                          WHERE taxa_index IN
                                          (SELECT thermo_index
                                          FROM taxa_pairs
                                          WHERE is_pair = True)""").df().values)
    size_tm_fafsa = int(con.execute("""SELECT COUNT(DISTINCT taxa_index)
                                           FROM fafsa_taxa
                                           WHERE taxa_index IN
                                           (SELECT DISTINCT meso_index
                                           FROM taxa_pairs)""").df().values)
    size_tt_fafsa = int(con.execute("""SELECT COUNT(DISTINCT taxa_index)
                                           FROM fafsa_taxa
                                           WHERE taxa_index IN
                                           (SELECT DISTINCT thermo_index
                                           FROM taxa_pairs)""").df().values)
    size_tm_small_diff = size_tm_16s_pair - size_tm_fafsa
    size_tt_small_diff = size_tt_16s_pair - size_tt_fafsa

    size_t_fafsa = size_tm_fafsa+size_tt_fafsa
    perc_t_no_pair = np.round(100*(size_tt_no_pair+size_tm_no_pair)/size_t_l2t, 1)
    perc_t_fafsa = np.round(100*size_t_fafsa/size_t_l2t, 1)
    perc_t_small_diff = np.round(100*(size_tt_small_diff+size_tm_small_diff)/size_tp_l2t, 1)

    # Builds and saves taxa Sankey in data folder
    fig2 = go.Figure(data=[go.Sankey(
    arrangement = 'snap',
    node = dict(
      pad = 20,
      thickness = 20,
      line = dict(color = 'black', width = 0.5),
      label = [f'learn2therm {str(np.round(size_t_l2t/1000, 1))+"k"}', 'Thermophile', 'Mesophile',
               '16S pair', '16S pair', f'No pairs ({perc_t_no_pair}%)',
               f'< {min_ogt_diff} \N{DEGREE SIGN}C diff ({perc_t_small_diff}%)',
               f"""FAFSA ({str(np.round(size_t_fafsa/1000, 1))+"k"},
               {perc_t_fafsa}%)"""],
      color = [sank_purple, sank_red, sank_blue, sank_red, sank_blue, sank_grey, sank_grey,
               sank_purple]
    ),
    link = dict(
      source = [2, 4, 0, 0, 1, 1, 2, 3, 3, 4],
      target = [5, 6, 1, 2, 3, 5, 4, 6, 7, 7],
      value = [size_tm_16s_pair, size_tm_small_diff, size_t_thermo, size_t_meso,
               size_tt_16s_pair, size_tt_no_pair, size_tm_16s_pair, size_tt_small_diff,
               size_tt_fafsa, size_tm_fafsa],
      color = [sank_grey_t, sank_grey_t, sank_red_t, sank_blue_t, sank_red_t, sank_grey_t,
               sank_blue_t, sank_grey_t, sank_red_t, sank_blue_t]
      ))])

    fig2.update_layout(title_text="Taxa Representation", font_size=16, font_family = 'Arial',
                       font_color = 'black')
    fig2.write_image(os.path.join('../data', newdir, 'taxa_sankey.png'), engine = 'kaleido',
                     scale = 6, width = 1280, height = 640)

    # Parameters for protein pair Sankey
    size_pp_l2t = int(con.execute("""SELECT COUNT(prot_pair_index)
                                     FROM protein_pairs""").df().values)
    size_pp_fafsa = int(con.execute("""SELECT COUNT(prot_pair_index)
                                           FROM fafsa_protein_pairs""").df().values)
    size_pp_small_diff = size_pp_l2t - size_pp_fafsa

    perc_pp_fafsa = np.round(100*size_pp_fafsa/size_pp_l2t, 1)
    perc_pp_small_diff = np.round(100*size_pp_small_diff/size_pp_l2t, 1)

    # Builds and saves protein pair Sankey in data folder
    fig3 = go.Figure(data=[go.Sankey(
    arrangement = 'snap',
    node = dict(
      pad = 25,
      thickness = 20,
      line = dict(color = "black", width = 0.5),
      label = [f'learn2therm ({str(np.round(size_pp_l2t/1000000, 1))+"m"})',
               f'< {min_ogt_diff} \N{DEGREE SIGN}C diff ({perc_pp_small_diff}%)',
               f"""FAFSA ({str(np.round(size_pp_fafsa/1000000, 1))+"m"},
               {perc_pp_fafsa}%)"""],
      color = [sank_purple, sank_grey, sank_purple]
    ),
    link = dict(
      source = [0, 0],
      target = [1, 2],
      value = [size_pp_small_diff, size_pp_fafsa],
      color = [sank_grey_t, sank_purple_t]
      ))])

    fig3.update_layout(title_text='Protein Pairs', font_family = 'Arial', font_size=16)
    fig3.write_image(os.path.join('../data', newdir, 'protein_pair_sankey.png'),
                     engine = 'kaleido', scale = 6, width = 1280, height = 640)

    # Parameters for protein Sankey
    size_p_l2t = int(con.execute("""SELECT COUNT(protein_int_index)
                                    FROM proteins""").df().values)
    size_p_meso = int(con.execute("""SELECT COUNT(DISTINCT meso_protein_int_index)
                                     FROM protein_pairs""").df().values)
    size_p_thermo = int(con.execute("""SELECT COUNT(DISTINCT thermo_protein_int_index)
                                       FROM protein_pairs""").df().values)
    size_pm_no_pair = int(con.execute("""SELECT COUNT(protein_int_index)
                                         FROM proteins
                                         WHERE taxa_index IN
                                         (SELECT meso_index
                                         FROM taxa_pairs
                                         WHERE is_pair = False)""").df().values)
    size_pt_no_pair = int(con.execute("""SELECT COUNT(protein_int_index)
                                         FROM proteins WHERE taxa_index IN
                                         (SELECT thermo_index
                                         FROM taxa_pairs
                                         WHERE is_pair = False)""").df().values)
    size_pm_16s_pair = int(con.execute("""SELECT COUNT(protein_int_index)
                                          FROM proteins
                                          WHERE taxa_index IN
                                          (SELECT meso_index FROM taxa_pairs
                                          WHERE is_pair = True)""").df().values)
    size_pt_16s_pair = int(con.execute("""SELECT COUNT(protein_int_index)
                                          FROM proteins
                                          WHERE taxa_index IN
                                          (SELECT thermo_index
                                          FROM taxa_pairs
                                          WHERE is_pair = True)""").df().values)
    size_pm_fafsa = int(con.execute("""SELECT COUNT(protein_int_index)
                                           FROM fafsa_proteins
                                           WHERE taxa_index in
                                           (SELECT DISTINCT meso_index
                                           FROM taxa_pairs)""").df().values)
    size_pt_fafsa = int(con.execute("""SELECT COUNT(protein_int_index)
                                           FROM fafsa_proteins
                                           WHERE taxa_index IN
                                           (SELECT DISTINCT thermo_index
                                           FROM taxa_pairs)""").df().values)
    size_pm_small_diff = size_pm_16s_pair - size_pm_fafsa
    size_pt_small_diff = size_pt_16s_pair - size_pt_fafsa
    size_p_null = int(con.execute("""SELECT COUNT(protein_int_index)
                                     FROM proteins
                                     WHERE taxa_index IN
                                     (SELECT DISTINCT taxa_index
                                     FROM taxa
                                     WHERE ogt IS NULL)""").df().values)
    size_p_fafsa = size_pm_fafsa+size_pt_fafsa

    perc_p_no_pair = np.round(100*(size_pt_no_pair+size_pm_no_pair)/size_p_l2t, 1)
    perc_p_fafsa = np.round(100*size_p_fafsa/size_p_l2t, 1)
    perc_p_small_diff = np.round(100*(size_pt_small_diff+size_pm_small_diff)/size_pp_l2t, 1)
    perc_p_null = np.round(100*size_p_null/size_p_l2t, 1)

    # Builds and saves protein Sankey in data folder
    fig4 = go.Figure(data=[go.Sankey(
    arrangement = 'snap',
    node = dict(
      pad = 50,
      thickness = 20,
      line = dict(color = 'black', width = 0.5),
      label = [f'learn2therm ({str(np.round(size_p_l2t/1000000, 1))+"m"})', 'Thermophile',
               'Mesophile', '16S pair', '16S pair', f'No pairs ({perc_p_no_pair}%)',
               f'< {min_ogt_diff} \N{DEGREE SIGN}C diff ({perc_p_small_diff}%)',
               f'FAFSA ({str(np.round(size_p_fafsa/1000000, 1))+"m"}, {perc_p_fafsa}%)',
               f'Null ({perc_p_null}%)'],
      color = [sank_purple, sank_red, sank_blue, sank_red, sank_blue, sank_grey, sank_grey,
               sank_purple, sank_grey]
    ),
    link = dict(
      source = [0, 0, 0, 1, 1, 2, 4, 4, 3, 3, 2],
      target = [8, 1, 2, 3, 5, 5, 6, 7, 6, 7, 4],
      value = [size_p_null, size_p_thermo, size_p_meso, size_pt_16s_pair, size_pt_no_pair,
               size_pm_no_pair, size_pm_small_diff, size_pm_fafsa, size_pt_small_diff,
               size_pt_fafsa, size_pm_16s_pair],
      color = [sank_grey_t, sank_red_t, sank_blue_t, sank_red_t, sank_grey_t, sank_grey_t,
               sank_grey_t, sank_blue_t, sank_grey_t, sank_red_t, sank_blue_t]
      ))])

    fig4.update_layout(title_text="Protein Representation", font_family = 'Arial', font_size=16)
    fig4.write_image(os.path.join('../data', newdir, 'protein_sankey.png'), engine = 'kaleido',
                     scale = 6, width = 1280, height = 640)
