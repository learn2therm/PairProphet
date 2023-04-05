#imports
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import os
import scipy.stats
import sklearn.preprocessing
import sklearn.model_selection
import sklearn.neighbors

df = pd.read_csv("/Users/loganroberts/Learn2Therm/ValidProt/data/learn2therm_sample_50k.csv")


#columns to drop. Determined during data exploration in Jupyter Notebook.
df = df.drop(columns = ['Unnamed: 0','thermo_index', 'm_protein_seq', 't_protein_seq', 'm_protein_desc', 't_protein_desc',
                        'query_align_cov_16s', 'subject_align_cov_16s', 'meso_index', 'meso_protein_int_index', 'local_gap_compressed_percent_id_16s', 
                        'scaled_local_query_percent_id_16s', 'scaled_local_symmetric_percent_id_16s','bit_score_16s'])


input_features = [columns for columns in df]
input_features.remove('bit_score')

target_feature = 'bit_score'


def split_data(dataframe):
    """
    Takes dataframe and splits it into dev and test sets.
    
    Params
    ----------
    dataframe: Pandas dataframe 

    Returns
    -------
    -Dev and test data
    
    """
    
    #test input data type
    if "pandas.core.frame.DataFrame" not in str(type(dataframe)):
        raise ValueError("Wrong input type!")
    else:
        pass
    
    #split data
    dev, test = sklearn.model_selection.train_test_split(dataframe, test_size=0.15, random_state=1)
    
    return dev, test


def train_reg(dev, test, columns = [],  target = []):
    """
    Takes dev and test dataframes and trains a standard Linear Regression model with selected data.
    
    Params
    ----------
    dev, test: Pandas dataframe previously split
    columns: list of strings, representing input features
    target: list of strings, representing target feature(s)

    Returns
    -------
    -Pearson correlation between each input and the output feature
    -Linear regression model
    -input test data vector (numpy array)
    -feature test vector (numpy array)
    """
    
    #test input arguments
    assert "pandas.core.frame.DataFrame" in str(type(dev))
    assert "pandas.core.frame.DataFrame" in str(type(test))
    assert "str" in str(type(columns[0]))
    assert "str" in str(type(target[0]))
    assert columns[0] in dev
    assert target in test
   
    #split into input and output feature(s)
    dev_X = dev[columns].values
    test_X = test[columns].values

    dev_y = dev[target].values.reshape(-1,1)
    test_y = test[target].values.reshape(-1,1)
    
    #return pearson correlation for linear correlation
    pearson_corr = dev.corr(method = 'pearson')[target]
    
    #scale data
    scaler = sklearn.preprocessing.StandardScaler()
    dev_X = scaler.fit_transform(dev_X)
    test_X = scaler.fit_transform(test_X)
    
    #train model
    model = sklearn.linear_model.LinearRegression()
    model = model.fit(dev_X, dev_y)
    
    return pearson_corr, model, test_X, test_y

p_corr, model, test_X, test_y = train_reg(split_data(df)[0], split_data(df)[1], columns = input_features, target = target_feature)


def test_reg(model, test_X, test_y):
    
    """
    Takes a trained model and test data and tests the model.
    
    Params
    ----------
    model: sklearn.linear_model
    test_X: numpy array
    test_y: numpy array

    Returns
    -------
    -Vector of predictions based on the model (numpy array)
    -prints R2 score of the model
    -prints mean absolute error and mean squared error from model
    """
    
    #test input arguments
    assert "sklearn.linear_model" in str(type(model))
    assert "numpy.ndarray" in str(type(test_X))
    assert "numpy.ndarray" in str(type(test_y))
    
    R2 = model.score(test_X, test_y)
    
    print('R2 score is: {}'.format(R2))

    preds = model.predict(test_X)

    mae = sklearn.metrics.mean_absolute_error(test_y, preds)
    mse = sklearn.metrics.mean_squared_error(test_y, preds)

    print("The MAE is : {} and the MSE is: {}".format(round(mae,6), round(mse,6)))
    
    return preds


def plot_regression(model, test_X, test_y):
    """
    Takes a test Linear Regression ML model and plots the predictions against actual values.
    
    Params
    ----------
    model: sklearn.linear_model
    test_X: numpy array
    test_y: numpy array

    Returns
    -------
    -Parity plot of predictions vs. observations
    -Equation of line of best fit (always linear)
    -R2 score
    """
    
    #test input arguments
    assert "sklearn.linear_model" in str(type(model))
    assert "numpy.ndarray" in str(type(test_X))
    assert "numpy.ndarray" in str(type(test_y))
    
    R2 = model.score(test_X, test_y)
    preds = test_reg(model, test_X, test_y)
    
    # make a plot to show the fit
    fig, ax = plt.subplots()

    # plot the actual vs. predicted values
    ax.scatter(test_y, preds, c='b',s=6, label='Parity Plot')
    plt.xlabel("Actual Values", fontsize=12)
    plt.ylabel("Predicted Values", fontsize=12)
    plt.title('Parity Plot', fontsize=15)

    #plot line of best fit
    slope, intercept = np.polyfit(np.ravel(test_y), preds, 1)
    plt.plot(np.ravel(test_y), slope*np.ravel(test_y) + intercept, c='green', lw=1)
    print('y = {}x + {}'.format(round(float(slope),4), round(float(intercept),4)))
    print('R2 = {}'.format(round(R2,4)))

