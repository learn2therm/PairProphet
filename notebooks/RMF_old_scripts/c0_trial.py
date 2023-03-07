from l2therm_in import path, sample_size, E_cutoff, ID_cutoff
import c1_preprocessing as c1

con = c1.connect_db(path)
df, time = c1.fetch_data(sample_size=sample_size, con = con, E_cutoff = E_cutoff, ID_cutoff = ID_cutoff)
print(df)
con.close()
