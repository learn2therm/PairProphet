import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import BaggingClassifier, \
    RandomForestClassifier, \
    AdaBoostClassifier, \
    GradientBoostingClassifier, \
    ExtraTreesClassifier

from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

names = ['LR', 'KNN', 'DT', 'NB', 'RF', 'Bagging', 'AB', 'GB', 'SVM']

#list of classifiers (hyperparameters optimized)
classifiers = [
    #Regression
    LogisticRegression(),
    #KNN (neighbors optimized iteratively)
    KNeighborsClassifier(n_neighbors=20),
    #Decision Tree
    DecisionTreeClassifier(max_features=None),
    #Gaussian
    GaussianNB(),
    #RF Classifier (with optuna)
    RandomForestClassifier(
        n_estimators=200,
        max_depth=None,
        max_samples=0.3,
        max_features=0.5,
        min_weight_fraction_leaf=0,
        min_samples_split=17),
    #RF Classifier with bagging (with optuna)
    BaggingClassifier(sklearn.ensemble.RandomForestClassifier
    (n_estimators=200, max_depth=None, 
     min_weight_fraction_leaf=0.000215),max_samples=0.5, 
    max_features=0.5),
    #AdaBoost (with optuna)
    AdaBoostClassifier(n_estimators=53, learning_rate=0.156),
    #Gradient Boosting (with optuna)
    GradientBoostingClassifier(n_estimators=100, learning_rate=1.0, 
                                max_depth=1),   
    #C-support vector classification (9)
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

def runClassifiers(filename:str, dataframe, columns=[], target=[], model=RandomForestClassifier()):

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
    -Accuracy score
    -area under ROC curve
    -train data (target)
    -validation data (features)
    -validation data (target)
    """
    
    dev, test = sklearn.model_selection.train_test_split(dataframe, test_size=0.15, random_state=1)

    train, val = sklearn.model_selection.train_test_split(dev, test_size=0.15, random_state=1)
    
    dev_X = dev[columns].values
    test_X = test[columns].values

    dev_y = dev[target].values.reshape(-1,1)
    test_y = test[target].values.reshape(-1,1) 

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

    Results = {}  # compare algorithms

    from sklearn.metrics import accuracy_score, \
        confusion_matrix, \
        roc_auc_score,\
        average_precision_score,\
        auc,\
        roc_curve, f1_score, recall_score, matthews_corrcoef, auc

    for classifier, name in zip(classifiers, names):
        accuracy = []
        avg_precision = []
        F1_Score = []
        AUC = []
        MCC = []
        Recall = []
        
        mean_TPR = 0.0
        mean_FPR = np.linspace(0, 1, 100)

        print('{} is done.'.format(classifier.__class__.__name__))

        model = classifier

        # model
        model.fit(train_X, train_y)

        preds = model.predict(val_X)

        # Calculate ROC Curve and Area the Curve
        proba_y = model.predict_proba(val_X)[:,1]
        FPR, TPR, _ = roc_curve(val_y, proba_y, pos_label=1)
        roc_auc = auc(FPR, TPR)
        
        #calculate scoring metrics
        #include option to return these scores
        accuracy.append(accuracy_score(y_pred=preds, y_true=val_y))
        avg_precision.append(average_precision_score(y_true=val_y, y_score=proba_y, pos_label=1))
        F1_Score.append(f1_score(y_true=val_y, y_pred=preds, pos_label=1))
        MCC.append(matthews_corrcoef(y_true=val_y, y_pred=preds))
        Recall.append(recall_score(y_true=val_y, y_pred=preds, pos_label=1))
        AUC.append(roc_auc)

        confusion_matrix = sklearn.metrics.confusion_matrix(y_pred=preds, y_true=val_y)
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
    
    F.close()
    
    return Results, model


# if __name__ == '__main__':
#     # print('Please, enter number of cross validation:')
#     import argparse
#     p = argparse.ArgumentParser(description='Run Machine Learning Classifiers.')

#     p.add_argument('-cv', '--nFCV', type=int, help='Number of crossValidation', default=10)
#     p.add_argument('-data', '--dataset', type=str, help='~/dataset.csv', default='optimumDataset.csv')
#     p.add_argument('-roc', '--auROC', type=int, help='Print ROC Curve', default=1, choices=[0, 1])
#     p.add_argument('-box', '--boxPlot', type=int, help='Print Accuracy Box Plaot', default=1, choices=[0, 1])

#     args = p.parse_args()

#     runClassifiers(args)