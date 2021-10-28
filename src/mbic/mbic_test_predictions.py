# AntiBiofilm Peptide Research
# Department of Computer Science and Engineering, Santa Clara University
# Author: Taylor Downey

# If SVM model decides peptide <=64uM then the SVR model is used to predict MBIC
# Hyperparameters for both models have already been tuned on the training set using cross-validation
# ------------------------------------------------------------------------------
#                               Libraries
# ------------------------------------------------------------------------------
import pandas as pd
import numpy as np
from sklearn import preprocessing
from sklearn.svm import SVR
from sklearn.svm import SVC
import json
from sklearn.decomposition import PCA

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
test_filename = '../../data/test_peptide_data.csv'
svm_features_filename = 'mbic_svm_forward_selection_features.json'
svr_features_filename = 'mbic_svr_forward_selection_features.json'
pred_filename = './mbic_predictions.csv'
split = 0.760   # Where to break peptides into class 0 and class 1

# Hyperparameters
svm_num_feats = 9
svm_pca_comp = 6
svm_c = 10
svm_g = 1000

svr_num_feats = 9
svr_pca_comp = 8
svr_c = 45
svr_g = 40

# ------------------------------------------------------------------------------
#                               Main
# ------------------------------------------------------------------------------
def main():
    
    # Prepare training peptides for SVM
    with open(svm_features_filename) as f:
        svm_feat_dict = json.load(f)
        svm_feat_dict = svm_feat_dict[0:svm_num_feats]
        
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

    # Prepare training peptides for SVM model
    min_max_scaler_svm = preprocessing.MinMaxScaler()
    X_norm_svm = min_max_scaler_svm.fit_transform(peptides_svm)
    pca_svm = PCA(n_components=svm_pca_comp)
    X_trans_svm = pca_svm.fit_transform(X_norm_svm) 
    SVC_rbf = SVC(kernel='rbf', C=svm_c, gamma=svm_g)
    print('Training SVM Peptides Shape: ', peptides_svm.shape)
    svm_fit = SVC_rbf.fit(X_trans_svm, y_svm)   # Train SVM model 

    # Prepare peptides for SVR
    with open(svr_features_filename) as f:
        svr_feat_dict = json.load(f)
        svr_feat_dict = svr_feat_dict[0:svr_num_feats]
        
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
    print('Training SVR Peptides Shape: ', peptides_svr.shape)
    svr_fit = SVR_rbf.fit(X_trans_svr, y_svr)   # Train SVR model 


    # Test peptides for SVM
    test_peptides_svm = pd.read_csv(test_filename)


    # Filter out columns based on popular feature list
    labels = test_peptides_svm.columns.values.tolist()
    for l in labels:
        if l not in svm_feat_dict:
            test_peptides_svm = test_peptides_svm.drop(columns=[l])

    # Prepare test peptides
    print('Test SVM Peptides Shape: ', test_peptides_svm.shape)
    X_norm_test_svm = min_max_scaler_svm.transform(test_peptides_svm)
    X_trans_tp_svm = pca_svm.transform(X_norm_test_svm)

    # Test peptides for SVR
    test_peptides_svr = pd.read_csv(test_filename)
    names = test_peptides_svr['Name'].tolist()
    dec_fuc = test_peptides_svr['Decision Fn'].tolist()

    # Filter out columns based on popular feature list
    labels = test_peptides_svr.columns.values.tolist()
    for l in labels:
        if l not in svr_feat_dict:
            test_peptides_svr = test_peptides_svr.drop(columns=[l])

    # Prepare test peptides
    print('Test SVR Peptides Shape: ', test_peptides_svr.shape)
    X_norm_test_svr = min_max_scaler_svr.transform(test_peptides_svr)
    X_trans_tp_svr = pca_svr.transform(X_norm_test_svr)

    # Train SVM model and then predict which bucket test peptides fall into
    test_peptide_classes = svm_fit.predict(X_trans_tp_svm)

    bucket0 = []
    for i,c in enumerate(test_peptide_classes):
        if c==1:
            bucket0.append(i)

    X_trans_tp_svr_small = X_trans_tp_svr[bucket0]
    test_peptide_mbic = svr_fit.predict(X_trans_tp_svr_small)

    mbic_names = []
    dec_funcs = []
    for j in bucket0:
        mbic_names.append(names[j])
        dec_funcs.append(dec_fuc[j])

    # Save MBIC predictions
    pred_results = list(zip(mbic_names, dec_funcs, test_peptide_mbic))
    df_pred_results = pd.DataFrame(pred_results, columns = ['Names', 'Decision Fn', 'Predicted MBIC Value']) 
    df_pred_results.to_csv(pred_filename, sep=',',index=False) 

if __name__ == "__main__":
    main()