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

from sklearn.metrics import accuracy_score, \
    confusion_matrix, \
    roc_auc_score,\
    average_precision_score,\
    auc,\
    roc_curve, f1_score, recall_score, matthews_corrcoef, auc

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

F = open('evaluationResults.txt', 'w')

F.write('Evaluation Scale:'+'\n')
F.write('0.0% <=Accuracy<= 100.0%'+'\n')
F.write('0.0 <=AUC<= 1.0'+'\n') #area under curve
F.write('0.0 <=auPR<= 1.0'+'\n')  # average_Precision
F.write('0.0 <=F1_Score<= 1.0'+'\n')
F.write('-1.0 <=MCC<= 1.0'+'\n')
F.write('_______________________________________'+'\n')


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

    accuracy = []
    avg_precision = []
    F1_Score = []
    AUC = []
    MCC = []
    Recall = []

    # test input arguments
    assert "sklearn" in str(type(model))
    assert "numpy.ndarray" in str(type(test_X))
    assert "numpy.ndarray" in str(type(test_y))

    # scale data
    scaler = sklearn.preprocessing.StandardScaler()
    test_X = scaler.fit_transform(test_X)

    #vector of predictions
    preds = model.predict(test_X)

    # calculate precision score
    precision_score = sklearn.metrics.precision_score(test_y, preds)

    # Calculate ROC Curve and Area the Curve
    proba_y = model.predict_proba(test_X)[:,1]
    FPR, TPR, _ = roc_curve(test_y, proba_y, pos_label=1)
    roc_auc = auc(FPR, TPR)
    
    #calculate scoring metrics
    #include option to return these scores
    accuracy.append(accuracy_score(y_pred=preds, y_true=test_y))
    avg_precision.append(average_precision_score(y_true=test_y, y_score=proba_y, pos_label=1))
    F1_Score.append(f1_score(y_true=test_y, y_pred=preds, pos_label=1))
    MCC.append(matthews_corrcoef(y_true=test_y, y_pred=preds))
    Recall.append(recall_score(y_true=test_y, y_pred=preds, pos_label=1))
    AUC.append(roc_auc)

    confusion_matrix = sklearn.metrics.confusion_matrix(y_pred=preds, y_true=test_y)
    sklearn.metrics.ConfusionMatrixDisplay(confusion_matrix).plot()

    accuracy = [_*100.0 for _ in accuracy]
    Results[name + ' Accuracy, F1 Score'] = [accuracy, F1_Score]
    
    F.write('Classifier: {}\n'.format(name))
    F.write('Accuracy: {0:.4f}%\n'.format(np.mean(accuracy)))
    F.write('AUC: {0:.4f}\n'.format( np.mean(AUC)))
    F.write('auPR: {0:.4f}\n'.format(np.mean(avg_precision))) # average_Precision
    F.write('F1_Score: {0:.4f}\n'.format(np.mean(F1_Score)))
    F.write('MCC: {0:.4f}\n'.format(np.mean(MCC)))

#         TN, FP, FN, TP = CM.ravel()
    F.write('Recall: {0:.4f}\n'.format( np.mean(Recall)) )
    F.write('_______________________________________'+'\n')

    return preds, precision_score


if __name__ == '__main__':
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
