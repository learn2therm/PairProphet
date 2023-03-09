"""
This module takes in a pandas dataframe from
c5_input_cleaning and runs it through a
RandomForestClassifier model from scitkit learn.
Returns a Boolean prediction for protein pair
functionality.
"""

import matplotlib.pyplot as plt
import pandas as pd
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

df['protein_match'] = df['t_protein_desc'].eq(df['m_protein_desc'])

# columns we don't want from our own database
# this will be cleaned in c5_input_cleaning.py
df = df.drop(
    columns=[
        'Unnamed: 0',
        'thermo_index',
        'm_protein_seq',
        't_protein_seq',
        'm_protein_desc',
        't_protein_desc',
        'query_align_cov_16s',
        'subject_align_cov_16s',
        'meso_index',
        'meso_protein_int_index',
        'local_gap_compressed_percent_id_16s',
        'scaled_local_query_percent_id_16s',
        'scaled_local_symmetric_percent_id_16s',
        'bit_score_16s',
        'm_ogt',
        't_ogt',
        'taxa_pair_index',
        'thermo_protein_int_index',
        'prot_pair_index',
        'ogt_difference'])


def train_model(dataframe, columns=[], target=[]):
    """
    Takes dataframe and splits it into a training and testing set.
    Note: Data is called train and test, but this test set is currently
    closer to a validation set. Keeping nomenclature to keep model robust.
    Trains a KNN classifier model with selected data.

    Params
    ----------
    dataframe: Pandas dataframe
    columns: list of strings, representing input features
    target: list of strings, representing target feature(s)

    Returns
    -------
    -Sk-learn model object
    -train data (features)
    -train data (target)
    -validation data (features)
    -validation data (target)
    """
    # split data
    dev, test = sklearn.model_selection.train_test_split(
        dataframe, test_size=0.15, random_state=1)

    # test input arguments
    assert "pandas.core.frame.DataFrame" in str(type(dev))
    assert "pandas.core.frame.DataFrame" in str(type(test))
    assert "str" in str(type(columns[0]))
    assert "str" in str(type(target[0]))

    # split into input and output feature(s)
    dev_X = dev[columns].values
    test_X = test[columns].values

    dev_y = dev[target].values.reshape(-1, 1)
    test_y = test[target].values.reshape(-1, 1)

    # scale data
    scaler = sklearn.preprocessing.StandardScaler()
    dev_X = scaler.fit_transform(dev_X)
    test_X = scaler.fit_transform(test_X)

    # train model
    model = sklearn.ensemble.RandomForestClassifier()
    model = model.fit(dev_X, dev_y.ravel())

    return model, dev_X, dev_y, test_X, test_y

# maybe combine with plotting model


def evaluate_model(model, test_X, test_y):
    """
    Takes a trained model and test data and tests the model.

    Params
    ----------
    model: sklearn.neighbors.KNeighborsClassifier
    test_X: numpy array
    test_y: numpy array

    Returns
    -------
    Vector of predictions based on the model (numpy array)
    """

    # test input arguments
    assert "sklearn" in str(type(model))
    assert "numpy.ndarray" in str(type(test_X))
    assert "numpy.ndarray" in str(type(test_y))

    preds = model.predict(test_X)

    return preds


def plot_model(model, test_X, test_y):
    """
    Takes a test classifier model and plots the confusion matrix.

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

    # test input arguments
    assert "sklearn" in str(type(model))
    assert "numpy.ndarray" in str(type(test_X))
    assert "numpy.ndarray" in str(type(test_y))

    score = model.score(test_X, test_y)
    preds = model.predict(test_X)

    # plot confusion matrix
    confusion_matrix = sklearn.metrics.confusion_matrix(preds, test_y)
    cm_plot = sklearn.metrics.ConfusionMatrixDisplay(confusion_matrix)

    cm_plot.plot(cmap=plt.cm.Blues)
    cm_plot.ax_.set_title('Confusion Matrix')

    return score


def rf_wrapper(dataframe):
    """
    Takes dataframe and runs it through kNN model.

    Params
    ----------
    dataframe: Pandas dataframe

    Returns
    -------
    -Target feature predictions
    -Parity plot
    """

    assert 'pandas.core.frame.DataFrame' in str(type(dataframe))

    # user inputs target feature
    target = 'protein_match'

    # define input features
    input_features = [columns for columns in dataframe]

    input_features.remove(target)

    # train the model based off data split
    model, dev_X, dev_y, test_X, test_y = train_model(
        dataframe, columns=input_features,
        target=target
    )

    # test the model and return predictions
    preds = evaluate_model(model, test_X, test_y)

    # plot the results of the model
    plot_model(model, test_X, test_y)

    return preds
