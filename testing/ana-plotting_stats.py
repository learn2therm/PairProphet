"""
TODO: Write a script that plots the statistics of the data
"""
# system dependencies

# library dependencies
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# local dependencies


####################
### PATHS & VARS ###
####################
# Paths

####################
#### SEABORN #######
####################
# seaborn settings
sns.set_context('talk')
sns.set_style('whitegrid')


####################
#### matplotlib ####
####################
# # matplotlib settings
# plt.rcParams.update({
#     'font.family': 'Arial',  # Times New Roman, Calibri
#     'font.weight': 'normal',
#     'mathtext.fontset': 'cm',
#     'font.size': 16,
#     'lines.linewidth': 2,
#     'axes.linewidth': 2,
#     'axes.spines.top': False,
#     'axes.spines.right': False,
#     'axes.titleweight': 'bold',
#     'axes.titlesize': 18,
#     'axes.labelweight': 'bold',
#     'xtick.major.size': 8,
#     'xtick.major.width': 2,
#     'ytick.major.size': 8,
#     'ytick.major.width': 2,
#     'figure.dpi': 80,
#     'legend.framealpha': 1, 
#     'legend.edgecolor': 'black',
#     'legend.fancybox': False,
#     'legend.fontsize': 14
# })

# ####################


####################
###    SCRIPT    ###
####################


# Read the data
## first number in variable name is e_value
## second number in variable name is Jaccard similarity score
data = pd.read_csv('./data/analysis/statistics_results.csv')



if __name__ == '__main__':
    # print(data.head())
    
    ###########
    # Summary #
    ###########
    print(data.describe())

    ###############
    # Manuplation #
    ###############

    # Define a parameter lambda
    lambda_param = 0.5  # This can be adjusted based on how much you want to penalize false_proportion

    # Compute the score
    data['score'] = data['true_proportion'] - lambda_param * data['false_proportion']
    print(f'new data head: {data.head()}')


    # ###########
    # # Heatmap #
    # ###########
    # # Compute correlation matrix
    # correlation_matrix_all = data[['e-value', 'jaccard_threshold', 'mean_acc_length', 'true_proportion', 'false_proportion']].corr()

    # # Plot heatmap
    # plt.figure(figsize=(10, 8))
    # sns.heatmap(correlation_matrix_all, annot=True, cmap='coolwarm', center=0)
    # plt.title('Correlation Heatmap of All Variables')
    # plt.savefig('./data/analysis/heatmap.png')  # Save the figure before showing it
    # plt.show()

    # ############
    # # Pairplot #
    # ############
    # # Create a pairplot to visualize the relationships between all pairs of variables
    # sns.pairplot(data, diag_kind='kde', plot_kws={'alpha': 0.6})
    # plt.savefig('./data/analysis/pairplot.png')  # Save the figure before showing it
    # plt.show()


    ################
    # Countour Plot #
    ################
    # Create a contour plot to visualize the relationship between e-value and jaccard_threshold

    # Define the space for interpolation
    grid_x, grid_y = np.mgrid[0:1:100j, data['e-value'].min():data['e-value'].max():100j]

    # Interpolate the data for score
    grid_z = griddata((data['jaccard_threshold'], data['e-value']), data['score'], (grid_x, grid_y), method='cubic')

    # Create contour plot
    plt.figure(figsize=(10, 8))
    cp = plt.contourf(grid_x, grid_y, grid_z, cmap='viridis')
    plt.colorbar(cp)
    plt.title('Score Contour Plot')
    plt.xlabel('Jaccard Threshold')
    plt.ylabel('e-value')
    plt.savefig('./data/analysis/cplot.png')  # Save the figure before showing it
    plt.show()

