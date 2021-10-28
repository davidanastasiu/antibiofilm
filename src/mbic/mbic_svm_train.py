# AntiBiofilm Peptide Research
# Department of Computer Science and Engineering, Santa Clara University
# Author: Taylor Downey

# A python script that performs forward selection and hyperparameter optimization
# in order to find the best performing SVM model for the MBIC peptides
# Loss function used is Matthew's Correlation Coefficient (MCC)
# The script dumps the features found from forward selection along with the 
# average MCC calculated with cross validation
# ------------------------------------------------------------------------------
#                               Libraries
# ------------------------------------------------------------------------------
import numpy as np
import pandas as pd
import json
import sys
import copy
import warnings
from sklearn import preprocessing
from sklearn.decomposition import PCA
from sklearn.utils.validation import column_or_1d
from sklearn.svm import SVC
from sklearn.metrics import matthews_corrcoef
from sklearn.model_selection import RepeatedStratifiedKFold
warnings.filterwarnings("ignore")

# ------------------------------------------------------------------------------
#                               Variables
# ------------------------------------------------------------------------------
training_filename = '../../data/mbic_training_data.csv'
svm_filename = 'mbic_svm_fs_mcc.txt'                       
svm_fs_features_filename = 'mbic_svm_forward_selection_features.json'

num_features = 200
C = [0.001, 0.01, 0.1, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 100, 1000]
gamma = [0.001, 0.01, 0.1, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 100]

# ------------------------------------------------------------------------------
#                               Main
# ------------------------------------------------------------------------------
def main():

    feat_list = []

    # Forward selection loop
    for nf in range(len(feat_list)+1, num_features + 1):

        peptides = pd.read_csv(training_filename)

        # Label peptides based on mbic values
        peptides.loc[(peptides['MBIC'] > 64), 'MBIC'] = 0  
        peptides.loc[(peptides['MBIC'] != 0), 'MBIC'] = 1

        y = peptides['MBIC'].to_numpy()

        peptides = peptides.drop(columns=['MBIC'])
        peptides = peptides.drop(columns=['Name'])
        peptides = peptides.drop(columns=['Pathogen'])
        peptides = peptides.drop(columns=['Type'])
        peptides = peptides.drop(columns=['Seq'])
        labels = peptides.columns.values.tolist()

        feat_MCC = []
        top_feat_list = []

        # Loop to populate list of features every iteration
        for l_num, l in enumerate(labels):
            small_feat_list = copy.deepcopy(feat_list)
            if l in feat_list:
                continue
            small_feat_list.append(l)
            peptides_small = peptides.copy()
            for f in labels:
                if f not in small_feat_list:
                    peptides_small = peptides_small.drop(columns=[f])

            num_feat = peptides_small.shape

            # Normalize features
            min_max_scaler = preprocessing.MinMaxScaler()
            X_norm = min_max_scaler.fit_transform(peptides_small)

            pca_features = []

            for j in range(1, num_feat[1]+ 1):
                pca_features.append(j)

            avg_MCC = []

            # Loops through each hyperparameter combination 
            for n in pca_features:
                pca = PCA(n_components=n)
                for c in C:
                    for g in gamma:
                        sys.stdout.write('Feature Loop: %s Feat num: %s PCA comp: %s c: %s \r' % (nf,l_num,n,c))
                        sys.stdout.flush()
                        X_trans = pca.fit_transform(X_norm)    
                        SVC_rbf = SVC(kernel='rbf', C=c, gamma=g)
                        MCC = []
                        rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats = 20)

                        # Cross validation
                        for train_index, test_index in rskf.split(X_trans, y):
                            X_train, X_test = X_trans[train_index], X_trans[test_index]
                            y_train, y_test = y[train_index], y[test_index]
                            y_train = y_train.reshape(-1,1)
                            y_train = column_or_1d(y_train, warn=False)
                            rbf_fit = SVC_rbf.fit(X_train, y_train)

                            y_pred = rbf_fit.predict(X_test)

                            mcc = matthews_corrcoef(y_test, y_pred)
                            MCC.append(mcc)

                        mcc_avg = np.mean(MCC)
                        avg_MCC.append(mcc_avg)

            # Take maximum MCC across all runs, select this feature, append to found list
            feat_MCC.append(np.max(avg_MCC))
            top_feat_list.append(copy.deepcopy(small_feat_list))
        
        idx = np.argmax(feat_MCC)
        feat_list = copy.deepcopy(top_feat_list[idx])

        # Dump MCC and current features found from forward selection
        with open (svm_filename, 'a', encoding="utf-8") as f:
                    f.write(str(nf) + '\t' + str(np.around(feat_MCC[idx], 3)) + '\n')

        with open(svm_fs_features_filename, 'w') as f:
            json.dump(feat_list, f) 

if __name__ == "__main__":
    main()

                        

       