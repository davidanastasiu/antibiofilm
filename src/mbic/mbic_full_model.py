# AntiBiofilm Peptide Research
# Department of Computer Science and Engineering, Santa Clara University
# Author: Taylor Downey

# A python script that uses the optimized hyperparameters found for both 
# the SVM and the SVR to create a prediction model
# Script prints the average RMSE of the full model when run with cross validation
# 
# NOTE: Given the small number of training samples available, the average RMSE 
# outputted will vary by about +- 5

# ------------------------------------------------------------------------------
#                               Libraries
# ------------------------------------------------------------------------------
import numpy as np
import pandas as pd
import json
import warnings
from sklearn import preprocessing
from sklearn.decomposition import PCA
from sklearn.utils.validation import column_or_1d
from sklearn.svm import SVC
from sklearn.svm import SVR
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import RepeatedStratifiedKFold
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
svm_features_filename = 'mbic_svm_forward_selection_features.json'
svr_features_filename = 'mbic_svr_forward_selection_features.json'
svr_svm_results = 'full_model_results.txt'

# Optimized Hyperparameters
svm_c = 10
svm_g = 1000
svm_pca_comp = 6
svm_num_feat = 9

svr_c = 45
svr_g = 40
svr_pca_comp = 8
svr_num_feat = 9

# ------------------------------------------------------------------------------
#                               Main
# ------------------------------------------------------------------------------
def main():

    # Prepare peptides for SVM
    with open(svm_features_filename) as f:
        svm_feat_dict = json.load(f)
        svm_feat_dict = svm_feat_dict[0:svm_num_feat]
        
    peptides_svm = pd.read_csv(training_filename)
    peptides_svm.loc[(peptides_svm['MBIC'] > 64), 'MBIC'] = 0  
    peptides_svm.loc[(peptides_svm['MBIC'] != 0), 'MBIC'] = 1

    # Filter out columns based on feat list
    labels = peptides_svm.columns.values.tolist()
    for l in labels:
        if l == 'MBIC':
            continue
        if l not in svm_feat_dict:
            peptides_svm = peptides_svm.drop(columns=[l])

    y_svm = peptides_svm['MBIC'].to_numpy()
    peptides_svm = peptides_svm.drop(columns=['MBIC'])

    min_max_scaler = preprocessing.MinMaxScaler()
    X_norm_svm = min_max_scaler.fit_transform(peptides_svm)
    pca_svm = PCA(n_components=svm_pca_comp)
    X_trans_svm = pca_svm.fit_transform(X_norm_svm) 
    SVC_rbf = SVC(kernel='rbf', C=svm_c, gamma=svm_g)

    # Prepare peptides for SVR
    with open(svr_features_filename) as f:
        svr_feat_dict = json.load(f)
        svr_feat_dict = svr_feat_dict[0:svr_num_feat]
        
    peptides_svr = pd.read_csv(training_filename)
    peptides_svr, _ = seperatePeptides(peptides_svr, 64)
        
    # Filter out columns based on feat list
    labels = peptides_svr.columns.values.tolist()
    for l in labels:
        if l == 'MBIC':
            continue
        if l not in svr_feat_dict:
            peptides_svr = peptides_svr.drop(columns=[l])

    y_svr = peptides_svr['MBIC'].to_numpy()
    peptides_svr = peptides_svr.drop(columns=['MBIC'])

    min_max_scaler_svr = preprocessing.MinMaxScaler()
    X_norm_svr = min_max_scaler_svr.fit_transform(peptides_svr)
    pca_svr = PCA(n_components=svr_pca_comp)
    X_trans_svr = pca_svr.fit_transform(X_norm_svr) 
    SVR_rbf = SVR(kernel='rbf', C=svr_c, gamma=svr_g)

    # Prepare test set of petides used by svr after training
    peptides_test_svr = pd.read_csv(training_filename)

    # Filter out columns based on feat list
    labels = peptides_test_svr.columns.values.tolist()
    for l in labels:
        if l == 'MBIC':
            continue
        if l not in svr_feat_dict:
            peptides_test_svr = peptides_test_svr.drop(columns=[l])

    y_svr2 = peptides_test_svr['MBIC'].to_numpy()
    peptides_test_svr = peptides_test_svr.drop(columns=['MBIC'])

    # Apply svr transformations on test set of peptides for svr
    X_norm_test_svr = min_max_scaler_svr.transform(peptides_test_svr)
    X_trans_test_svr = pca_svr.transform(X_norm_test_svr)
        
    # Cross validation applied to full model
    rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats = 20)
    RMSE = []
    cnt = 1
    for train_index, test_index in rskf.split(X_trans_svm, y_svm):  
        X_train, X_test = X_trans_svm[train_index], X_trans_svm[test_index]
        y_train, y_test = y_svm[train_index], y_svm[test_index]
        y_train = y_train.reshape(-1,1)
        y_train = column_or_1d(y_train, warn=False)
        svm_fit = SVC_rbf.fit(X_train, y_train)
        y_pred = svm_fit.predict(X_test)
        
        train_index_svr = []
        test_index_svr = []
        y_train_svr = []
        y_test_svr = []
        for i in range(0, len(y_train)):
            if(y_train[i] == 0):
                continue
            else:
                train_index_svr.append(train_index[i])
                
        X_train_svr = X_trans_svr[train_index_svr]
        y_train_svr = y_svr[train_index_svr]
        svr_fit = SVR_rbf.fit(X_train_svr, y_train_svr)
        
        y_train_svr = []
        for i in range(0, len(y_pred)):
            if(y_pred[i] == 0):
                continue
            else:
                test_index_svr.append(test_index[i])
        
        X_test_svr = X_trans_test_svr[test_index_svr]
        y_test_svr = y_svr2[test_index_svr]
        y_pred_svr = SVR_rbf.predict(X_test_svr)
        
        rmse = np.sqrt(mean_squared_error(y_test_svr, y_pred_svr))

        cnt = cnt + 1
        
        with open (svr_svm_results, 'a', encoding="utf-8") as sfile:
            sfile.write(str(rmse) + '\n')
        
        RMSE.append(rmse)
            
    rmse_avg = np.average(RMSE)
            
    print('RMSE average: ' + str(rmse_avg))

if __name__ == "__main__":
    main()
    
    
