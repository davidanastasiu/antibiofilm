# AntiBiofilm Peptide Research
# Department of Computer Science and Engineering, Santa Clara University
# Author: Taylor Downey

# A python script that performs forward selection and hyperparameter optimization
# in order to find the best performing SVR model for the MBIC peptides
# Loss function used is Root Mean Squared Error (RMSE)
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
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import RepeatedKFold
warnings.filterwarnings("ignore")

# ------------------------------------------------------------------------------
#                               Functions
# ------------------------------------------------------------------------------
def seperatePeptides(peptides, threshold):

    columns = ['MBIC']
    
    filterMBIC = (peptides[columns] <= threshold).all(axis=1)  
    lower_peptides = peptides[filterMBIC]

    filterMBIC = (peptides[columns] > threshold).all(axis=1)  
    upper_peptides = peptides[filterMBIC]
    
    return lower_peptides, upper_peptides

# ------------------------------------------------------------------------------
#                               Variables
# ------------------------------------------------------------------------------
training_filename = '../../data/mbic_training_data.csv'
svr_filename = 'mbic_svr_fs_rmse.txt'                       
svr_fs_features_filename = 'mbic_svr_forward_selection_features.json'

num_features = 200
C = [0.001, 0.01, 0.1, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 100, 1000]
gamma = [0.001, 0.01, 0.1, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 100]

# ------------------------------------------------------------------------------
#                               Main
# ------------------------------------------------------------------------------
def main():

    feat_list = []

    # Forward selection loop
    for nf in range(len(feat_list) + 1, num_features + 1):

        peptides = pd.read_csv(training_filename)

        # Get lower valued peptides
        peptides, _ = seperatePeptides(peptides, 64)

        y = peptides['MBIC'].to_numpy()

        peptides = peptides.drop(columns=['MBIC'])
        peptides = peptides.drop(columns=['Name'])
        peptides = peptides.drop(columns=['Pathogen'])
        peptides = peptides.drop(columns=['Type'])
        peptides = peptides.drop(columns=['Seq'])
        labels = peptides.columns.values.tolist()

        feat_RMSE = []
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

            avg_RMSE = []

            # Loops through each hyperparameter combination
            for n in pca_features:
                pca = PCA(n_components=n)
                for c in C:
                    for g in gamma:
                        sys.stdout.write('Feature Loop: %s Feat num: %s PCA comp: %s c: %s \r' % (nf,l_num,n,c))
                        sys.stdout.flush()
                        X_trans = pca.fit_transform(X_norm)    
                        SVR_rbf = SVR(kernel='rbf', C=c, gamma=g)
                        RMSE = []
                        rskf = RepeatedKFold(n_splits=5, n_repeats = 20)

                        # Cross validation
                        for train_index, test_index in rskf.split(X_trans, y):
                            X_train, X_test = X_trans[train_index], X_trans[test_index]
                            y_train, y_test = y[train_index], y[test_index]
                            y_train = y_train.reshape(-1,1)
                            y_train = column_or_1d(y_train, warn=False)
                            rbf_fit = SVR_rbf.fit(X_train, y_train)

                            y_pred = rbf_fit.predict(X_test)
                    
                            rmse = np.sqrt(mean_squared_error(y_test, y_pred))
                            RMSE.append(rmse)

                        rmse_avg = np.mean(RMSE)
                        avg_RMSE.append(rmse_avg)

            # Take minimum RMSE across all runs, select this feature, append to found list
            feat_RMSE.append(np.min(avg_RMSE))
            top_feat_list.append(copy.deepcopy(small_feat_list))
           
        idx = np.argmin(feat_RMSE)
        feat_list = copy.deepcopy(top_feat_list[idx])

        # Dump RMSE and current features found from forward selection
        with open (svr_filename, 'a', encoding="utf-8") as f:
                    f.write(str(nf) + '\t' + str(np.around(feat_RMSE[idx], 3)) + '\n')

        with open(svr_fs_features_filename, 'w') as f:
            json.dump(feat_list, f) 

if __name__ == "__main__":
    main()