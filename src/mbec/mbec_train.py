# AntiBiofilm Peptide Research
# Department of Computer Science and Engineering, Santa Clara University
# Author: Taylor Downey

# A python script that performs forward selection and hyperparameter optimization
# in order to find the best performing SVR model for the MBEC peptides
# Loss function used is RMSE
# The script dumps the features found from forward selection along with the 
# average RMSE calculated with cross validation
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
from sklearn.svm import SVR
from sklearn.feature_selection import RFE
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import RepeatedKFold
warnings.filterwarnings("ignore")

# ------------------------------------------------------------------------------
#                                Variables
# ------------------------------------------------------------------------------
num_features = 200
C = [0.001, 0.01, 0.1, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 100, 1000]
gamma = [0.001, 0.01, 0.1, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 100]
# ------------------------------------------------------------------------------
#                                   Main
# ------------------------------------------------------------------------------
def main():

    feat_list = []

    # Forward selection loop
    for n in range(0, num_features + 1):
        print('\nFeatures: ' + str(n))
        peptides = pd.read_csv('../../data/mbec_training_data.csv')
        peptides = peptides.drop(columns=['Name', 'Seq'])
        labels = peptides.columns.values.tolist()
        labels.remove('MBEC(uM)')
        feat_RMSE = []
        top_feat_list = []

        # Loop to populate list of features every iteration
        for i, l in enumerate(labels):
            small_feat_list = copy.deepcopy(feat_list)
            if l in feat_list:
                continue
            small_feat_list.append(l)
            peptides_small = peptides.copy()
            for f in labels:
                if f not in small_feat_list:
                    peptides_small = peptides_small.drop(columns=[f])

            y = peptides_small['MBEC(uM)'].to_numpy()
            peptides_small = peptides_small.drop(columns=['MBEC(uM)'])
            num_feat = peptides_small.shape

            # Normalize features
            min_max_scaler = preprocessing.MinMaxScaler()
            X_norm = min_max_scaler.fit_transform(peptides_small)

            pca_features = []

            for j in range(1, num_feat[1]+ 1):
                pca_features.append(j)
            
            avg_RMSE = []

            # Loops through each hyperparameter combination 
            for n in pca_features:
                for c in C:
                    sys.stdout.write('Feature Loop: %s PCA comp: %s c: %s \r' % (i,n,c))
                    sys.stdout.flush()
                    for g in gamma:
                        SVR_rbf = SVR(kernel='rbf', C=c, gamma=g)
                        RMSE = []

                        kf = RepeatedKFold(n_splits=5, n_repeats = 20)

                        # Cross validation
                        for train_index, test_index in kf.split(X_norm):
                            X_train, X_test = X_norm[train_index], X_norm[test_index]
                            y_train, y_test = y[train_index], y[test_index]
                            y_train = y_train.reshape(-1,1)
                            y_train = column_or_1d(y_train, warn=False)

                            model = SVR_rbf.fit(X_train, y_train)
                            y_pred = model.predict(X_test)

                            test_error = mean_squared_error(y_test, y_pred)
                            RMSE.append(np.sqrt(test_error))

                        avg_RMSE.append(np.mean(RMSE))

            # Take minimum RMSE across all runs, select this feature, append to found list
            feat_RMSE.append(np.min(avg_RMSE))
            top_feat_list.append(copy.deepcopy(small_feat_list))

        # Dump RMSE and current features found from forward selection
        idx = np.argmin(feat_RMSE)
        feat_list = copy.deepcopy(top_feat_list[idx])
        with open ('mbec_fs_rmse.txt', 'a', encoding="utf-8") as f:
                    f.write(str(n) + '\t' + str(np.around(feat_RMSE[idx], 3)) + '\n')

        with open('mbec_forward_selection_features.json', 'w') as f:
            json.dump(feat_list, f)          

if __name__ == "__main__":
    main()