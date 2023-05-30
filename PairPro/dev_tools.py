from PairPro.preprocessing import connect_db
import pandas as pd

def build_sample_l2t(db_in, db_out, size):
    '''
    Generates a sample l2t relational database of given size. Note that the final
    size will about 30% of 'size' due to pair filtering.

    Args:
        db_in (str): Path to full size l2t database
        db_out (int): Path to sample l2t database to be created
        size (str): Number of pairs to sample for test database. 
                    Final size will be about 30% of this.

    Returns:
        None. Database file is saved at db_out.

    Raises:
        None.
    '''
    con, _ = connect_db(db_in)    
    
    # Make sample protein pairs table
    cmd1 = f"""CREATE TEMP TABLE samp_protein_pairs AS
             SELECT * FROM protein_pairs 
             USING SAMPLE {size}"""

    con.execute(cmd1)
    
    # Make sample proteins table
    cmd2 = """CREATE TEMP TABLE samp_proteins AS
             SELECT * FROM proteins
             WHERE proteins.pid IN 
             (SELECT DISTINCT meso_pid FROM samp_protein_pairs)
             OR proteins.pid IN
             (SELECT DISTINCT thermo_pid FROM samp_protein_pairs)"""

    con.execute(cmd2)
    
    # Make sample taxa table
    cmd3 = """CREATE TEMP TABLE samp_taxa AS
             SELECT * FROM taxa
             WHERE taxa.taxid IN 
             (SELECT DISTINCT meso_taxid FROM samp_protein_pairs)
             OR taxa.taxid IN
             (SELECT DISTINCT thermo_taxid FROM samp_protein_pairs)"""

    con.execute(cmd3)
    
    # Make sample taxa_pairs table
    cmd4 = """CREATE TEMP TABLE samp_taxa_pairs AS
             SELECT * FROM taxa_pairs
             WHERE taxa_pairs.query_id IN 
             (SELECT DISTINCT taxid FROM samp_taxa)
             AND taxa_pairs.subject_id IN 
             (SELECT DISTINCT taxid FROM samp_taxa)"""

    con.execute(cmd4)
    
    # Make sample taxa_lab table
    cmd5 = """CREATE TEMP TABLE samp_taxa_pairs_lab AS
             SELECT * FROM taxa_pairs_lab
             WHERE taxa_pairs_lab.__index_level_0__ IN 
             (SELECT __index_level_0__ FROM samp_taxa_pairs)
          """

    con.execute(cmd5)
    
    # Grab new tables as df and close connection to large database
    samp_protein_pairs = con.execute("""SELECT * FROM samp_protein_pairs""").df()
    samp_proteins = con.execute("""SELECT * FROM samp_proteins""").df()
    samp_taxa = con.execute("""SELECT * FROM samp_taxa""").df()
    samp_taxa_pairs = con.execute("""SELECT * FROM samp_taxa_pairs""").df()
    samp_taxa_pairs_lab = con.execute("""SELECT * FROM samp_taxa_pairs_lab""").df()
    con.close()
    
    con2 = duckdb.connect(db_out)
    con2.execute("""CREATE OR REPLACE TABLE protein_pairs AS SELECT * FROM samp_protein_pairs""")
    con2.execute("""CREATE OR REPLACE TABLE proteins AS SELECT * FROM samp_proteins""")
    con2.execute("""CREATE OR REPLACE TABLE taxa AS SELECT * FROM samp_taxa""")
    con2.execute("""CREATE OR REPLACE TABLE taxa_pairs AS SELECT * FROM samp_taxa_pairs""")
    con2.execute("""CREATE OR REPLACE TABLE taxa_pairs_lab AS SELECT * FROM samp_taxa_pairs_lab""")
    
    con2.close()