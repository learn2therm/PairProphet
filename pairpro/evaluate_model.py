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


def evaluate_model(output_path, model, dataframe):
    '''
    Takes a trained model and test data and tests the model.

    Args:
        output path: File path: str
        model: sklearn.neighbors.KNeighborsClassifier
        test_X: numpy array
        test_y: numpy array

    Returns:
        Vector of predictions (numpy arrray)
    '''
    from sklearn.metrics import accuracy_score, \
    confusion_matrix, \
    roc_auc_score,\
    average_precision_score,\
    auc,\
    roc_curve, f1_score, recall_score, matthews_corrcoef, auc
    
    # initialize empty eval results file
    F = open('evaluationResults.txt', 'w')

    F.write('Evaluation Scale:' + '\n')
    F.write('0.0% <=Accuracy<= 100.0%' + '\n')
    F.write('0.0 <=AUC<= 1.0' + '\n')  # area under curve
    F.write('0.0 <=auPR<= 1.0' + '\n')  # average_Precision
    F.write('0.0 <=F1_Score<= 1.0' + '\n')
    F.write('-1.0 <=MCC<= 1.0' + '\n')
    F.write('_______________________________________' + '\n')

    model = model

    df_seqs = dataframe[['m_protein_seq', 't_protein_seq', 'prot_pair_index']]

    dataframe = dataframe.drop(columns=['m_protein_seq', 't_protein_seq',
                                        'prot_pair_index'])

    target = 'protein_match'
    features = [columns for columns in dataframe]
    features.remove(target)
    print(features)

    # split into input and output feature(s)
    test_X = dataframe[features].values
    test_y = dataframe[target].values.reshape(-1, 1)

    # scale data
    scaler = sklearn.preprocessing.StandardScaler()
    test_X = scaler.fit_transform(test_X)

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

    # vector of predictions
    preds = model.predict(test_X)

    # calculate precision score
    precision_score = sklearn.metrics.precision_score(test_y, preds)

    # Calculate ROC Curve and Area the Curve
    proba_y = model.predict_proba(test_X)[:, 1]
    FPR, TPR, _ = roc_curve(test_y, proba_y, pos_label=1)
    roc_auc = auc(FPR, TPR)

    # calculate scoring metrics
    # include option to return these scores
    accuracy.append(accuracy_score(y_pred=preds, y_true=test_y))
    avg_precision.append(
        average_precision_score(
            y_true=test_y,
            y_score=proba_y,
            pos_label=1))
    F1_Score.append(f1_score(y_true=test_y, y_pred=preds, pos_label=1))
    MCC.append(matthews_corrcoef(y_true=test_y, y_pred=preds))
    Recall.append(recall_score(y_true=test_y, y_pred=preds, pos_label=1))
    AUC.append(roc_auc)

    confusion_matrix = sklearn.metrics.confusion_matrix(
        y_pred=preds, y_true=test_y)
    sklearn.metrics.ConfusionMatrixDisplay(confusion_matrix).plot()

    accuracy = [_*100.0 for _ in accuracy]

    F.write('Accuracy: {0:.4f}%\n'.format(np.mean(accuracy)))
    F.write('AUC: {0:.4f}\n'.format(np.mean(AUC)))
    F.write(
        'auPR: {0:.4f}\n'.format(
            np.mean(avg_precision)))  # average_Precision
    F.write('F1_Score: {0:.4f}\n'.format(np.mean(F1_Score)))
    F.write('MCC: {0:.4f}\n'.format(np.mean(MCC)))

#         TN, FP, FN, TP = CM.ravel()
    F.write('Recall: {0:.4f}\n'.format(np.mean(Recall)))
    F.write('_______________________________________' + '\n')

    # merge dataframes together to report results
    df_seqs['prediction'] = preds

    # for multi-class classifier problem
    # preds_df = pd.DataFrame(preds, columns=['hmmer_pred', 'structure_pred'])
    # df_seqs['hmmer_pred'] = preds_df['hmmer_pred'].values
    # df_seqs['structure_pred'] = preds_df['structure_pred'].values
    # df_seqs['hmmer_strucutre_match'] = df_seqs['hmmer_pred'] == df_seqs['structure_pred']
    
    # merge dataframes together to report results (1 class)
    df_seqs['prediction'] = preds

    # save to csv
    df_seqs.to_csv(f'{output_path}predictions.csv')

    return preds, precision_score, df_seqs