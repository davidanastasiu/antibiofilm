# AntiBiofilm Peptide Research
# Department of Computer Science and Engineering, Santa Clara University
# Author: Taylor Downey

# A python script that builds a model with the features found during forward 
# selection and the optimized hyperparameters. 
# Model predicts the mbec value of each given test peptide
# ------------------------------------------------------------------------------
#                               Libraries
# ------------------------------------------------------------------------------
import pandas as pd
import numpy as np
import json
from sklearn import preprocessing
from sklearn.svm import SVR
from sklearn.decomposition import PCA

# ------------------------------------------------------------------------------
#                               Variables
# ------------------------------------------------------------------------------
training_filename = '../../data/mbec_training_data.csv'
test_filename = '../../data/test_peptide_data.csv'
fs_filename = 'forward_selection_features.json'    # Features chosen using Forward Selection
pred_filename = './mbec_predictions.csv'

# Optimized hyperparameters
num_feats = 12
pca_comp = 7
c = 1000
g = 20

# ------------------------------------------------------------------------------
#                               Main
# ------------------------------------------------------------------------------

def main():

    # Training peptides
    training_peptides = pd.read_csv(training_filename)
    s = training_peptides.shape

    # Load forward selection features
    with open(fs_filename) as f:
        feat_dict = json.load(f)

    feat_dict = feat_dict[0:num_feats]

    # Filter out all features not chosen by forward selection
    labels = training_peptides.columns.values.tolist()
    for l in labels:
        if l == 'MBEC(uM)':
            continue
        if l not in feat_dict:
            training_peptides = training_peptides.drop(columns=[l])

    # Test peptides
    test_peptides = pd.read_csv(test_filename)

    # Load forward selection features
    with open(fs_filename) as f:
        feat_dict = json.load(f)

    feat_dict = feat_dict[0:num_feats]

    names = test_peptides['Name'].tolist()
    dec_fuc = test_peptides['Decision Fn'].tolist()

    # Filter out features all features not chosed by forward selection
    labels = test_peptides.columns.values.tolist()
    for l in labels:
        if l not in feat_dict:
            test_peptides = test_peptides.drop(columns=[l])

    rows, cols = test_peptides.shape
    mbec_predictions = np.zeros((3, rows))  # numpy array to store mbec predictions

    # Prepare training peptides 
    min_max_scaler = preprocessing.MinMaxScaler()
    y = training_peptides['MBEC(uM)'].to_numpy()                        # Convert response variable to array
    training_peptides = training_peptides.drop(columns=['MBEC(uM)'])    # Drop MBEC from training features
    print('Training Peptides Shape: ', training_peptides.shape)
    X_norm = min_max_scaler.fit_transform(training_peptides)            # Scale training data

    pca = PCA(n_components=pca_comp)
    X_trans = pca.fit_transform(X_norm)     # Perform PCA on training data

    # Prepare test peptides
    print('Test Peptides Shape: ', test_peptides.shape)
    X_norm_test = min_max_scaler.transform(test_peptides)
    X_trans_tp = pca.transform(X_norm_test)

    # Build SVR model and train
    SVC_rbf = SVR(kernel='rbf', C=c, gamma=g)
    rbf_fit = SVC_rbf.fit(X_trans, y)

    # Predict test peptides
    y_test_pred = rbf_fit.predict(X_trans_tp)

    # Save MBEC predictions
    pred_results = list(zip(names, dec_fuc, y_test_pred))
    df_pred_results = pd.DataFrame(pred_results, columns = ['Names', 'Decision fn', 'Predicted MBEC']) 
    df_pred_results.to_csv(pred_filename, sep=',',index=False) 

if __name__ == "__main__":
    main()