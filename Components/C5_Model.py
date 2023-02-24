import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import os
import scipy.stats
import sklearn.preprocessing
import sklearn.model_selection
import sklearn.neighbors

df = pd.read_csv("/Users/loganroberts/Learn2Therm/ValidProt/data/Sample.csv")

df = df.drop(columns = ['Unnamed: 0', 'prot_pair_index', 'meso_seq', 'thermo_seq'])

def train_reg(dataframe, columns = [],  target = []):
    """
    Takes a dataframe and trains a standard Linear Regression model with selected data.
    
    Params
    ----------
    dataframe: Pandas dataframe
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
    assert "pandas.core.frame.DataFrame" in str(type(dataframe))
    assert "str" in str(type(columns[0]))
    assert "str" in str(type(target[0]))
    
    #split data
    dev, test = sklearn.model_selection.train_test_split(dataframe, test_size=0.20, random_state=1)
    
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
    

p_corr, model, test_X, test_y = train_reg(df, columns = ['meso_ogt', 'thermo_ogt', 'scaled_local_symmetric_percent_id',
                                                         'local_E_value','local_gap_compressed_percent_id'], target = ['scaled_local_query_percent_id'])

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
    
    
    preds = test_reg(model, test_X, test_y)
    R2 = model.score(test_X, test_y)
    
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
