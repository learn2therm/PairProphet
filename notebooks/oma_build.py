import duckdb
import sys
import time

sys.path.append('/home/ryfran/PairProphet/pairpro')

import utils

# WARNING: This code will take a long time to run and probably requires ~32 GB of memory. The files are huge.
# Expect resulting db file to be > 100 GB.
t1 = time.time()
DB_DIR = '/mnt/s/OMA/'
DB_NAME = 'oma.db'

pair_file = 'oma-pairs.txt'
seq_file = 'oma-seq.fa'
uniprot_file = 'oma-uniprot.txt'
prok = '/mnt/s/OMA/prokaryotes.cdna.fa'
parquet_dir = '/mnt/s/OMA/parquet'

# This code should install fasql ad hoc. It allows duckdb to read fasta files
con = duckdb.connect('/mnt/s/OMA/db/oma.db', config={'allow_unsigned_extensions': True})
con.execute("SET custom_extension_repository='dbe.wheretrue.com/fasql/latest';")
con.execute("INSTALL fasql;")
con.execute("LOAD fasql;")

# Make parquet files
pairs_cols = ['protein1', 'protein2', 'mult', 'group']
uniprot_cols = ['oma_id', 'uniprot_id']
utils.split_txt(in_path=f'{DB_DIR}{pair_file}', out_dir=f'{parquet_dir}/pairs', cols=pairs_cols)
utils.split_txt(in_path=f'{DB_DIR}{uniprot_file}', out_dir=f'{parquet_dir}/uniprot', cols=uniprot_cols)

# Make pairs SQL tables from parquets
utils.parqs_to_db(parq_dir=f'{parquet_dir}/pairs', table_name='pairs', con=con, cols=['protein1', 'protein2'])
utils.parqs_to_db(parq_dir=f'{parquet_dir}/uniprot', table_name='uniprot', con=con)
t2 = time.time()
print(f'Completed parquet build: {t2-t1} seconds')

con.execute(f"""CREATE OR REPLACE TABLE proteins AS SELECT description, sequence FROM read_fasta('{DB_DIR}{seq_file}')""")

t3 = time.time()
print(f'Completed proteins build: {t3-t2} seconds')

# Creates a new table with only prokaryotic protein pairs.
con.execute(f"""CREATE OR REPLACE TABLE prok_pairs AS SELECT
                protein1 AS protein1_oma_id,
                protein2 AS protein2_oma_id,
                FROM pairs 
                WHERE protein1 IN (SELECT description FROM read_fasta('{prok}'))
                AND protein2 IN (SELECT description FROM read_fasta('{prok}'))""")

# Updates the prok_pairs table with uniprot ids. This and the previous steps must be split to avoid memory issues.
con.execute("""CREATE OR REPLACE TABLE prok_pairs AS SELECT
               protein1_oma_id,
               protein2_oma_id,
               u1.uniprot_id AS protein1_uniprot_id,
               u2.uniprot_id AS protein2_uniprot_id,
               FROM prok_pairs
               LEFT JOIN uniprot u1 ON u1.oma_id = protein1_oma_id
               LEFT JOIN uniprot u2 ON u2.oma_id = protein2_oma_id""")

# At this pont, db file should be functional.
t4 = time.time()

print(f'Completed OMA build: {t4-t3} seconds')
con.close()