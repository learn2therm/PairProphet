import pandas as pd
import user_blast as ub

data = pd.read_csv('../notebooks/c0-c2_exploration_plotting_sampling/test100.csv')[['m_protein_seq','t_protein_seq']]
newdf = ub.make_blast_df(data)
print(newdf)