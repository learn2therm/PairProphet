import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import joblib


# testing function

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

def evaluate_model(output_path, model, target: list, dataframe):
    '''
    Takes a trained model and test data and tests the model.
    Runs a single or multi-class Classifier depending on input.

    Args:
        output path: File path: str
        model: sklearn.neighbors.KNeighborsClassifier
        target: target for classifier (list)
        dataframe: pandas dataframe

    Returns:
        Vector of predictions (numpy arrray)
        precision score (numpy array)
        results (csv)
    '''
    from sklearn.metrics import accuracy_score, \
    confusion_matrix, \
    roc_auc_score,\
    average_precision_score,\
    auc,\
    roc_curve, f1_score, recall_score, matthews_corrcoef, auc
    
    if 'structure_match' not in target:
        # initialize empty eval results file
        F = open('evaluationResults.txt', 'w')

        F.write('Evaluation Scale:' + '\n')
        F.write('0.0% <=Accuracy<= 100.0%' + '\n')
        F.write('0.0 <=AUC<= 1.0' + '\n')  # area under curve
        F.write('0.0 <=auPR<= 1.0' + '\n')  # average_Precision
        F.write('0.0 <=F1_Score<= 1.0' + '\n')
        F.write('-1.0 <=MCC<= 1.0' + '\n')
        F.write('_______________________________________' + '\n')

        results_df = dataframe[['m_protein_seq', 't_protein_seq']]

        dataframe = dataframe.drop(columns=['m_protein_seq', 't_protein_seq'])

        features = [columns for columns in dataframe.drop(columns=target)]

        # split into input and output feature(s)
        test_X = dataframe[features].values
        test_y = dataframe[target].values.reshape(-1, 1)

        # scale data
        scaler = sklearn.preprocessing.StandardScaler()
        test_X = scaler.fit_transform(test_X)

        # test input arguments
        assert "sklearn" in str(type(model))
        assert "numpy.ndarray" in str(type(test_X))
        assert "numpy.ndarray" in str(type(test_y))

        # vector of predictions
        preds = model.predict(test_X)

        # Calculate ROC Curve and Area the Curve
        proba_y = model.predict_proba(test_X)[:, 1]
        FPR, TPR, _ = roc_curve(test_y, proba_y, pos_label=1)
        roc_auc = auc(FPR, TPR)

        # calculate scoring metrics
        # include option to return these scores
        accuracy = 100*(accuracy_score(y_pred=preds, y_true=test_y))
        avg_precision = average_precision_score(
                y_true=test_y,
                y_score=proba_y,
                pos_label=1)
        F1_Score = f1_score(y_true=test_y, y_pred=preds, pos_label=1)
        MCC = matthews_corrcoef(y_true=test_y, y_pred=preds)
        Recall = recall_score(y_true=test_y, y_pred=preds, pos_label=1)
        AUC = roc_auc

        # confusion_matrix = sklearn.metrics.confusion_matrix(
        #     y_pred=preds, y_true=test_y)
        # sklearn.metrics.ConfusionMatrixDisplay(confusion_matrix).plot()

        F.write('Accuracy: {0:.4f}%\n'.format(accuracy))
        F.write('AUC: {0:.4f}\n'.format(AUC))
        F.write(
            'auPR: {0:.4f}\n'.format(avg_precision)
            )  # average_Precision
        F.write('F1_Score: {0:.4f}\n'.format(F1_Score))
        F.write('MCC: {0:.4f}\n'.format(MCC))

    #         TN, FP, FN, TP = CM.ravel()
        F.write('Recall: {0:.4f}\n'.format(Recall))
        F.write('_______________________________________' + '\n')

        # merge dataframes together to report results
        results_df['prediction'] = preds

        # save to csv
        results_df.to_csv('predictions.csv')

    else:
         # initialize empty eval results file
        F = open('evaluationResults.txt', 'w')

        F.write('Evaluation Scale:' + '\n')
        F.write('0.0% <=Accuracy<= 100.0%' + '\n')
        F.write('0.0 <=AUC<= 1.0' + '\n')  # area under curve
        F.write('0.0 <=auPR<= 1.0' + '\n')  # average_Precision
        F.write('0.0 <=F1_Score<= 1.0' + '\n')
        F.write('-1.0 <=MCC<= 1.0' + '\n')
        F.write('_______________________________________' + '\n')

        results_df = dataframe[['m_protein_seq', 't_protein_seq']]

        dataframe = dataframe.drop(columns=['m_protein_seq', 't_protein_seq'])

        features = [columns for columns in dataframe.drop(columns=target)]

        # split into input and output feature(s)
        test_X = dataframe[features].values
        test_y = dataframe[target].values

        # scale data
        scaler = sklearn.preprocessing.StandardScaler()
        test_X = scaler.fit_transform(test_X)

        # test input arguments
        assert "sklearn" in str(type(model))
        assert "numpy.ndarray" in str(type(test_X))
        assert "numpy.ndarray" in str(type(test_y))

        # vector of predictions
        preds = model.predict(test_X)
        hmmer_preds = preds[:, 0]
        structure_preds = preds[:, 1]
        assert len(hmmer_preds) == len(structure_preds)

        # calculate scoring metrics
        accuracy = 100*accuracy_score(y_pred=preds, y_true=test_y)
        Precision = precision_score(
                y_pred=preds,
                y_true=test_y,
                average=None,
                pos_label=1
                )
        F1_Score = f1_score(y_true=test_y, y_pred=preds, pos_label=1, average=None)
        Recall = recall_score(y_true=test_y, y_pred=preds, pos_label=1, average=None)

        print(accuracy)
        print(Precision)
        print(F1_Score)
        print(Recall)
    
        # confusion_matrix = sklearn.metrics.confusion_matrix(
        #     y_pred=preds, y_true=test_y)
        # sklearn.metrics.ConfusionMatrixDisplay(confusion_matrix).plot()

        F.write('Accuracy: {0:.4f}%\n'.format(accuracy))
        F.write('Mean Precision: {0:.4f}\n'.format(np.mean(Precision)))
        F.write('Mean F1_Score: {0:.4f}\n'.format(np.mean(F1_Score)))

    #         TN, FP, FN, TP = CM.ravel()
        F.write('Mean Recall: {0:.4f}\n'.format(np.mean(Recall)))
        F.write('_______________________________________' + '\n')

        # merge dataframes together to report results
        results_df['hmmer_prediction'] = hmmer_preds
        results_df['structure_prediction'] = structure_preds
        results_df['hmmer_structure_match'] = results_df['hmmer_prediction'] == results_df['structure_prediction']

        # save to csv
        results_df.to_csv('predictions.csv')

    return preds, results_df