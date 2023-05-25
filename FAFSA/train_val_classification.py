"""
This module takes in a pandas dataframe from
c5_input_cleaning and runs it through a
RandomForestClassifier model from scitkit learn.
Returns a Boolean prediction for protein pair
functionality.
"""

import matplotlib.pyplot as plt
import sklearn.preprocessing
import sklearn.model_selection
import sklearn.neighbors
import sklearn.ensemble
import sklearn.feature_selection
import sklearn.metrics


def train_model(dataframe, columns=[], target=[]):
    """
    Takes dataframe and splits it into a training and testing set. 
    Trains a RF Classifier with data.

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
    train, val = sklearn.model_selection.train_test_split(
        dataframe, test_size=0.15, random_state=1)

    # test input arguments
    assert "pandas.core.frame.DataFrame" in str(type(train))
    assert "pandas.core.frame.DataFrame" in str(type(val))
    assert "str" in str(type(columns[0]))
    assert "str" in str(type(target[0]))

    # split into input and output feature(s)
    train_X = train[columns].values
    val_X = val[columns].values

    train_y = train[target].values.reshape(-1, 1)
    val_y = val[target].values.reshape(-1, 1)

    # scale data
    scaler = sklearn.preprocessing.StandardScaler()
    train_X = scaler.fit_transform(train_X)
    val_X = scaler.fit_transform(val_X)

    # train model with hyperparams optimized
    model = sklearn.ensemble.RandomForestClassifier(
        n_estimators=200,
        max_depth=None,
        max_samples=0.3,
        max_features=0.5,
        min_weight_fraction_leaf=0,
        min_samples_split=17)

    model = model.fit(train_X, train_y.ravel())

    return model, train_X, train_y, val_X, val_y


def evaluate_model(model, val_X, val_y):
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
    assert "numpy.ndarray" in str(type(val_X))
    assert "numpy.ndarray" in str(type(val_y))

    preds = model.predict(val_X)

    # not printed during model validation step
    precision_score = sklearn.metrics.precision_score(val_y, preds)

    return preds, precision_score


def plot_model(model, val_X, val_y):
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
    assert "numpy.ndarray" in str(type(val_X))
    assert "numpy.ndarray" in str(type(val_y))

    score = model.score(val_X, val_y)
    preds = model.predict(val_X)

    # plot confusion matrix
    confusion_matrix = sklearn.metrics.confusion_matrix(preds, val_y)
    cm_plot = sklearn.metrics.ConfusionMatrixDisplay(confusion_matrix)

    cm_plot.plot(cmap=plt.cm.Blues)
    cm_plot.ax_.set_title('Confusion Matrix')

    return score


def rf_wrapper(dataframe):
    """
    Takes dataframe and runs it through RandomForest model.

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
    model, _, _, val_X, val_y = train_model(
        dataframe, columns=input_features,
        target=target
    )

    # test the model and return predictions
    preds, _ = evaluate_model(model, val_X, val_y)

    # plot the results of the model
    score = plot_model(model, val_X, val_y)

    return preds, score, model
