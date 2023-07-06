import duckdb
import pandas as pd
import sys
import time

# Test file for sampler function. 

sys.path.append('/home/ryfran/PairProphet/pairpro')

import utils

con = duckdb.connect('/mnt/s/OMA/db/oma.db')

for sample in [100, 10000, 1000000]:
    print(f'Starting sample {sample}')
    t1 = time.time()
    df = utils.oma_sample(con=con, size=sample)
    df.to_csv(f'/mnt/s/OMA/samples/oma_sample_{sample}.csv', index=False)
    t2 = time.time()
    print(f'Completed sample {sample}: {t2-t1} seconds')