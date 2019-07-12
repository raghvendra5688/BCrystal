import xgboost as xgb
import numpy as np
import pandas as pd
from sklearn.utils import shuffle
from sklearn.utils import class_weight


## Load the features
df = pd.read_csv('features.csv')
X_train = []
for i in range(len(df)):
    df.append(df.values[i][1:])

## Load the Labels
df = pd.read_csv('./Data/Train/Train_True_Labels.csv')
y_train = []
for i in range(len(df)):
    y_train.append(df.values[i][0])

X_train = np.array(X_train)
y_train = np.array(y_train)

parmas = {
'max_depth':5,
'min_child_weight':6,
'eta':0.1,
'subsample':0.7,
'colsample_bytree':0.6,
'objective':'binary:logistic'
}

## Assigning weights for each class because our data set is imbalanced.
weights = np.zeros(len(y_train))
weights[y_train1 == 0] = 1
weights[y_train1 == 1] = 2.94
dtrain = xgb.DMatrix(X_train, label=y_train, weight=weights)

bs = xgb.train(parmas, dtrain)

bs.save_model('train.model')




