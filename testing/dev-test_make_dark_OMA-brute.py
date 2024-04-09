import pandas as pd
import numpy as np
import time
import pairpro.user_blast as ub
from denseweight import DenseWeight

def bad_sampler(df, size, cols):
    bad_sample = pd.DataFrame(columns=cols)
    print(f'starting with size {size}')
    while bad_sample.shape[0] < size:
        df_it = pd.DataFrame([df.sample(frac=.1).reset_index(drop=True)['protein1_seq'], df.sample(frac=.1).reset_index(drop=True)['protein2_seq']]).T
        df_it.columns = ['query', 'subject']
        out, _ = ub.make_blast_df(df_it, path = None)
        bad_sample = pd.concat([bad_sample, out], ignore_index=True, axis=0)
        print(f'added {out.shape[0]} rows, total {bad_sample.shape[0]}')
    return bad_sample

def weighted_filter(x, weight, lowpass=True):
    new_x = []
    norm = max(weight(x))
    if lowpass == True:
        minloc = np.where(weight(x) == min(weight(x)))
        minval = min(weight(x))/norm


        for i, val in enumerate(x):

            if val > x[list(minloc[0])[0]]:
                if np.random.rand() < weight(val)/norm:
                    new_x.append(i)
            else:
                if np.random.rand() < minval:
                    new_x.append(i)

    else:
        for i, val in enumerate(x):
            if np.random.rand() < weight(val)/norm:
                    new_x.append(i)
            
    return new_x

size = 1000
df = pd.read_csv('/mnt/s/OMA/OMA_seq_1m_limit.csv')
df.drop(columns=['Unnamed: 0', 'protein1_uniprot_id', 'protein2_uniprot_id'], inplace=True)

scratch = df.sample(size).reset_index(drop=True)

t1 = time.time()
print('Pair blast')

if size > 500:
    pairs = df.sample(500).reset_index(drop=True)
    out, _ = ub.make_blast_df(pairs, path = None)
    print(f'Number of rows: {out.shape[0]}')

    while out.shape[0] < size:
        pairs = df.sample(500).reset_index(drop=True)
        out_it, _ = ub.make_blast_df(pairs, path = None)
        out = pd.concat([out, out_it.reset_index(drop=True)], ignore_index=True)
        print(f'Number of rows: {out.shape[0]}')

else:        
    pairs = df.sample(size).reset_index(drop=True)
    out, _ = ub.make_blast_df(pairs, path = None)

t2 = time.time()
print('Time: ', t2-t1)

print('Scratch blast')
nonpairs = bad_sampler(scratch, size*100, cols = out.columns)

t3 = time.time()

print('Time: ', t3-t2)

out['Pair'] = True

scratch_filtered = nonpairs[nonpairs['local_gap_compressed_percent_id'] > 0]
scratch_filtered['Pair'] = False

dw = DenseWeight(alpha=1)
dw.fit(list(out['local_gap_compressed_percent_id']))

idx = weighted_filter(scratch_filtered['local_gap_compressed_percent_id'], dw)
scratch_weighted = scratch_filtered.iloc[idx]

out = pd.concat([out, scratch_weighted], ignore_index=True)

out.to_csv('/mnt/s/OMA/delete.csv')