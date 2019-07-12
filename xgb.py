#Code to run for testing on new protein sequences

import matplotlib as mpl
import matplotlib.pyplot as plt
from xgboost import XGBClassifier
import xgboost as xgb
import numpy as np
import pandas as pd
np.random.seed(1337)
import io
import pickle
import csv
import argparse
import os
from Bio import SeqIO
import shap

parser = argparse.ArgumentParser()
parser.add_argument('file', type=argparse.FileType('r'))
parser.add_argument('file1', type=argparse.FileType('r'))
parser.add_argument('directory', help='the destination')
parser.add_argument("--encoding", default="utf_8")   
args = parser.parse_args()

csv_reader = csv.reader(
    args.file,
    delimiter=",",
    quotechar='"'
)

Seq_id = []
with args.file1 as file:
        fasta_sequences = SeqIO.parse(file,'fasta')
        for i,fasta in enumerate(fasta_sequences, 1):
            if (i <= 5):
                name, sequence = fasta.id, str(fasta.seq)
                Seq_id.append(name)
            else:
                break

#Get and Manipulate Sequences
x = []
for row in csv_reader:
    x.append(row)

X_test1 = []
for i in range(1,len(x)):
    X_test1.append(x[i][1:])

df = pd.read_csv('Feature_Names.csv')
X1 = []
for i in range(len(df.values)):
    X1.append(df.values[i][0])

X_test = []
X_test.append(X1)
for i in range(len(X_test1)):
    X_test.append(np.array(X_test1[i]).tolist())

#Add feature name information on test set
df = pd.DataFrame(X_test)
df.to_csv("Metadata_Test.csv", header=False)
df = pd.read_csv('Metadata_Test.csv')
feature_names = []
for i in range(1, len(df.columns)):
    feature_names.append(df.columns[i])
X_test2 = df[feature_names]

X_test1 = X_test2
X_test2 = np.array(X_test2)
dtest = xgb.DMatrix(X_test2)

#Run the model on new test set
bs = xgb.Booster({'nthread': 4})  # init model
bs.load_model('xgb_cross_best_score_15_less90.model')  # load data
yprob = bs.predict(dtest)

scores = []
for i in range(len(yprob)):
    scores.append(yprob[i])

out_directory = 'mkdir '+args.directory
if os.path.exists(args.directory) is False:
	os.system(out_directory)
	
#Output the results in <output-dir>/prediction.csv file
with open(args.directory+'/prediction.csv', 'w') as f:
    write = csv.writer(f)
    write.writerows(zip(Seq_id, yprob))
f.close()

#Plot the top 10 features driving crystallization propensity for each protein in test file
colors = ["#7fffd4","#ff0000"]
for i in range(0, (len(Seq_id))):
    data_for_prediction = X_test1.iloc[[i]]
    data_for_prediction.columns = feature_names
    explainer = shap.TreeExplainer(bs)
    shap_values = explainer.shap_values(data_for_prediction.values)
    data_for_prediction = np.round(data_for_prediction, 3)
    shap.initjs()
    fig = shap.force_plot(explainer.expected_value, shap_values, data_for_prediction, plot_cmap=colors, matplotlib=True, show=False)
    plt.figure(figsize=(11.5, 8), dpi=300)

    x2 = np.argsort(np.abs(shap_values[0]))
    x2 = x2[::-1]
    yaxis = []
    for j in range(10):
        yaxis.append(shap_values[0][x2[j]])
    
    xax2 = []
    for j in range(10):
        xax2.append(feature_names[x2[j]] + "=" + str(np.round(X_test2[i][x2[j]], 3)))
    yaxis = yaxis[::-1]
    xaxis = xax2[::-1]
    yaxis_pos = []
    yaxis_neg = []
    xaxis_pos = []
    xaxis_neg = []
    colors = []
    labels = []
    for j in range(len(yaxis)):
        if(yaxis[j] >=0):
            colors.append("#7fffd4")
            yaxis_pos.append(yaxis[j])
            xaxis_pos.append(xaxis[j])
            labels.append("-> Crys")
        else:
            colors.append("#ff0000")
            yaxis_neg.append(yaxis[j])
            xaxis_neg.append(xaxis[j])
            labels.append("-> Non-Crys")
    if len(yaxis_pos) > 0 & len(yaxis_neg) > 0  :
        plt.barh(xaxis,yaxis, color=colors, label=labels)
        plt.legend(["-> Crys","-> Non-Crys"],loc=4)
    if len(yaxis_neg) > 0 & len(yaxis_pos)==0:
        plt.barh(xaxis_neg,yaxis_neg, color="#ff0000",label="-> Non-Crys")
        plt.legend(loc=4)
    if len(yaxis_pos)> 0 & len(yaxis_neg)==0:
        plt.barh(xaxis_pos,yaxis_pos, color="#7fffd4",label="-> Crys")
        plt.legend(loc=4)
    plt.xlim(-1.3, 1)
    plt.xlabel("SHAP Values")
    plt.title('Top Features for: '+Seq_id[i]+', score = '+str(yprob[i]))
    plt.savefig(args.directory+'/bar_plot_'+str(i+1)+'.png', bbox_inches='tight')
	
os.system('rm Metadata_Test.csv')
print("DONE")
