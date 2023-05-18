import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# scikit-learn :
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import BaggingClassifier, \
    RandomForestClassifier, \
    AdaBoostClassifier, \
    GradientBoostingClassifier, \
    ExtraTreesClassifier

# create dictionary of models
names = ['LR', 'KNN', 'DT', 'NB', 'RF', 'Bagging', 'AB', 'GB', 'SVM']

# list of classifiers (hyperparameters optimized)
classifiers = [
    # Regression
    LogisticRegression(),
    # KNN (neighbors optimized iteratively)
    KNeighborsClassifier(n_neighbors=20),
    # Decision Tree
    DecisionTreeClassifier(max_features=None),
    # Gaussian
    GaussianNB(),
    # RF Classifier (with optuna)
    RandomForestClassifier(
        n_estimators=200,
        max_depth=None,
        max_samples=0.3,
        max_features=0.5,
        min_weight_fraction_leaf=0,
        min_samples_split=17),
    # RF Classifier with bagging (with optuna)
    BaggingClassifier(RandomForestClassifier
                      (n_estimators=200, max_depth=None,
                       min_weight_fraction_leaf=0.000215), max_samples=0.5,
                      max_features=0.5),
    # AdaBoost (with optuna)
    AdaBoostClassifier(n_estimators=53, learning_rate=0.156),
    # Gradient Boosting (with optuna)
    GradientBoostingClassifier(n_estimators=100, learning_rate=1.0,
                               max_depth=1),
    # C-support vector classification (9)
    #     SVC(),
]


def test_model(args, test_X, test_y):
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

    model = args.classifier

    # test input arguments
    assert "sklearn" in str(type(model))
    assert "numpy.ndarray" in str(type(test_X))
    assert "numpy.ndarray" in str(type(test_y))

    preds = model.predict(test_X)

    # not printed during model validation step
    precision_score = sklearn.metrics.precision_score(test_y, preds)

    return preds, precision_score

if __name__ == '__main__':
     # print('Please, enter number of cross validation:')
    import argparse
    p = argparse.ArgumentParser(description='Choose your classification model!')

    p.add_argument('-classifier', '--classifier', type=str, help='Specify which classifier you want to use.', default=RandomForestClassifier(
        n_estimators=200,
        max_depth=None,
        max_samples=0.3,
        max_features=0.5,
        min_weight_fraction_leaf=0,
        min_samples_split=17))
    
    args = p.parse_args()

    test_model(args, test_X, test_y)
