import pandas as pd
import numpy as np

import c1

path = '/mnt/c/Users/Ryan/Downloads/learn2therm_sample_50k.zip'

df = c1.fetch_data(path = path, form = 'csv')

print(df.head())