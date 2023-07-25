"""
TODO: Write a script that plots the statistics of the data
"""
# system dependencies

# library dependencies
import matplotlib.pyplot as plt
import pandas as pd

# local dependencies


####################
### PATHS & VARS ###
####################
# Paths


####################
#### matplotlib ####
####################
# matplotlib settings
plt.rcParams.update({
    'font.family': 'Helvetica',  # Times New Roman, Calibri
    'font.weight': 'normal',
    'mathtext.fontset': 'cm',
    'font.size': 16,
    'lines.linewidth': 2,
    'axes.linewidth': 2,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.titleweight': 'bold',
    'axes.titlesize': 18,
    'axes.labelweight': 'bold',
    'xtick.major.size': 8,
    'xtick.major.width': 2,
    'ytick.major.size': 8,
    'ytick.major.width': 2,
    'figure.dpi': 80,
    'legend.framealpha': 1, 
    'legend.edgecolor': 'black',
    'legend.fancybox': False,
    'legend.fontsize': 14
})

####################


####################
###    SCRIPT    ###
####################


# Read the data
## first number in variable name is e_value
## second number in variable name is Jaccard similarity score
data_1_1 = pd.read_csv('./data/analysis/1_statistics_results.csv')
# data_1_2 = pd.read_csv('./data/analysis/2_statistics_results.csv')
# data_1_3 = pd.read_csv('./data/analysis/3_statistics_results.csv')
# data_1_4 = pd.read_csv('./data/analysis/4_statistics_results.csv')
# data_2_1 = pd.read_csv('./data/analysis/5_statistics_results.csv')
# data_2_2 = pd.read_csv('./data/analysis/6_statistics_results.csv')
# data_2_3 = pd.read_csv('./data/analysis/7_statistics_results.csv')
# data_2_4 = pd.read_csv('./data/analysis/8_statistics_results.csv')
# data_3_1 = pd.read_csv('./data/analysis/9_statistics_results.csv')
# data_3_2 = pd.read_csv('./data/analysis/10_statistics_results.csv')
# data_3_3 = pd.read_csv('./data/analysis/11_statistics_results.csv')
# data_3_4 = pd.read_csv('./data/analysis/12_statistics_results.csv')
# data_4_1 = pd.read_csv('./data/analysis/13_statistics_results.csv')
# data_4_2 = pd.read_csv('./data/analysis/14_statistics_results.csv')
# data_4_3 = pd.read_csv('./data/analysis/15_statistics_results.csv')
# data_4_4 = pd.read_csv('./data/analysis/16_statistics_results.csv')


if __name__ == '__main__':
    print(data_1_1.head())