import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.stats
import sklearn.preprocessing
import sklearn.model_selection
import sklearn.neighbors
import sklearn.ensemble
import sklearn.feature_selection


# learn how to repeat this
df_original = pd.read_csv('learn2therm_sample_50k.csv')

"""
This cleaning step is done because we want to use thermophyllic protein
description as our target. I am eliminating categories without significant
value counts in order to build a better model with our "dummy" target. 
This code will be deleted when we get our real target from Component 3.
"""

categories = df_original['t_protein_desc'].value_counts()
categories = categories.iloc[categories.values > 500]
categories_dict = {item: None for item in categories.index}
list_of_cats = list(categories_dict.keys())
list_of_cats

# process dataframe
df = df_original[df_original.t_protein_desc.isin(list_of_cats)]


#columns we don't want from our own database
df = df.drop(columns=['Unnamed: 0', 'thermo_index', 'm_protein_seq',
                't_protein_seq', 'm_protein_desc', 'query_align_cov_16s',
                'subject_align_cov_16s', 'meso_index', 'meso_protein_int_index',
                'local_gap_compressed_percent_id_16s', 'scaled_local_query_percent_id_16s',
                'scaled_local_symmetric_percent_id_16s', 'bit_score_16s', 'm_ogt', 
                't_ogt', 'taxa_pair_index', 'thermo_protein_int_index', 
                'prot_pair_index', 'ogt_difference'])


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


def train_reg(dev, test, columns=[],  target=[]):
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
   
    #split into input and output feature(s)
    dev_X = dev[columns].values
    test_X = test[columns].values

    dev_y = dev[target].values.reshape(-1, 1)
    test_y = test[target].values.reshape(-1, 1)
    
    #scale data
    scaler = sklearn.preprocessing.StandardScaler()
    dev_X = scaler.fit_transform(dev_X)
    test_X = scaler.fit_transform(test_X)
    
    #train model
    model = sklearn.neighbors.KNeighborsClassifier(n_neighbors=20)
    model = model.fit(dev_X, dev_y.ravel())
    
    return model, dev_X, dev_y, test_X, test_y


def test_reg(model, test_X, test_y):
    
    """
    Takes a trained model and test data and tests the model.
    
    Params
    ----------
    model: sklearn.neighbors.KNeighborsClassifier
    test_X: numpy array
    test_y: numpy array

    Returns
    -------
    -Vector of predictions based on the model (numpy array)
    -prints R2 score of the model
    -prints mean absolute error and mean squared error from model
    -prints mean, standard deviation, and relative standard deviation for the target feature
    """
    
    #test input arguments
    assert "sklearn" in str(type(model))
    assert "numpy.ndarray" in str(type(test_X))
    assert "numpy.ndarray" in str(type(test_y))
    
    score = model.score(test_X, test_y)
    
    print('Model score is: {}'.format(round(score, 4)))
    
    preds = model.predict(test_X)
    
    return preds


def plot_regression(model, test_X, test_y):
    """
    Takes a test KNN Classifier model and plots the confusion matrix.
    
    Params
    ----------
    model: sklearn.neighbors.KNeighborsClassifier
    test_X: numpy array
    test_y: numpy array

    Returns
    -------
    -Confusion predictions vs. observations
    -Model score
    """
    
    #test input arguments
    assert "sklearn" in str(type(model))
    assert "numpy.ndarray" in str(type(test_X))
    assert "numpy.ndarray" in str(type(test_y))
    
    score = model.score(test_X, test_y)
    preds = test_reg(model, test_X, test_y)
   
    # plot confusion matrix
    confusion_matrix = sklearn.metrics.confusion_matrix(preds, test_y)
    cm_plot = sklearn.metrics.ConfusionMatrixDisplay(confusion_matrix)
    
    cm_plot.plot(cmap=plt.cm.Blues)
    cm_plot.ax_.set_title('Confusion Matrix')


def kNN_wrapper(dataframe):
    """
    Takes dataframe and runs it through kNN model.
    
    Params
    ----------
    dataframe: Pandas dataframe 

    Returns
    -------
    -Pearson correlation
    -Target feature predictions
    -Parity plot
    
    """
    
    assert 'pandas.core.frame.DataFrame' in str(type(dataframe))

    #user inputs target feature
    target = 't_protein_desc'
    
    #define input features
    input_features = [columns for columns in dataframe]
    
    input_features.remove(target)
    
    print(input_features)

    #train the model based off data split
    model, dev_X, dev_y, test_X, test_y = train_reg(
                                        split_data(dataframe)[0], split_data(dataframe)[1], 
                                        columns=input_features, target=target
                                                   )
    
    #test the model and return predictions
    preds = test_reg(model, test_X, test_y)
    
    #plot the results of the model
    plot_regression(model, test_X, test_y)
    
    return preds

#define function for Jensen-Shannon Divergence 
def JSD_dev_and_test(data_sample_1, data_sample_2, axis=None):
    """Jensen-Shannon divergence between two data samples.
    
    Formatted for use with scipy bootstrapping. Creates a vector of values covering
    the data range, the fits scipy historgram PDFs to both data samples. Each PDF is
    evaluated at the vector of values, and those PDF values are used to compute
    JSD. Note that becaue the empirical distribution is a histogram, JSD will rarely
    be less than about ~0.1 due to discretization artifacts.
    
    Parameters
    ----------
    data_sample_1 : numpy array or pandas series
        observations in sample 1
    data_sample_2 : numpy array or pandas series
        observations in sample 2
    axis : ignored
    
    Returns
    -------
    float : Value of JSD
    """
    # get the min and max values to compute PDF
    global_minimum = min([data_sample_1.min(), data_sample_2.min()])
    global_maximum = min([data_sample_1.max(), data_sample_2.max()])
    width = global_maximum - global_minimum
    
    x = np.linspace(global_minimum - 0.1*width, global_maximum + 0.1*width, 1000)
    
    # create empirical PDFs for the samples
    s1_distro = scipy.stats.rv_histogram(np.histogram(data_sample_1))
    s2_distro = scipy.stats.rv_histogram(np.histogram(data_sample_2))
    
    # compute JSD
    JSD = scipy.spatial.distance.jensenshannon(s1_distro.pdf(x), s2_distro.pdf(x))
    return [JSD]