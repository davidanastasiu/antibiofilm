# Identification of distinct characterisitcs of antibiofilm peptides and prospection of diverse sources for efficacious sequences

This respository contains data and code associated with the paper found [here](https://biorxiv.org/cgi/content/short/2021.09.28.462235v1) 

## Authors 

Bipasa Bose, Taylor Downey, Anand K. Ramasubramanian, and David C. Anastasiu

### Summary

In this work, we developed machine learning models to identify the distinguishing characteristics
of known antibiofilm peptides, and to mine peptide databases from diverse habitats to classify new peptides
with potential antibiofilm activities. Additionally, we used the reported minimum inhibitory/eradication
concentration (MBIC/MBEC) of the antibiofilm peptides to create a regression model on top of the classification
model to predict the effectiveness of new antibiofilm peptides. We used a positive dataset containing 242
antibiofilm peptides, and a negative dataset which, unlike previous datasets, contains peptides that are 
likely to promote biofilm formation. We utilized our classification-regression pipeline to evaluate 135,015 
peptides from diverse sources and identified antibiofilm peptide candidates that are efficacious against 
preformed biofilms at micromolar concentration.

### Data

MBIC training data can be found [here](/data/mbic_training_data.csv)  
MBEC training data can be found [here](/data/mbec_training_data.csv)  
Test peptides can be found [here](/data/test_peptide_data.csv)  

### MBIC/MBEC scripts

Two independent models were developed in order to predict the minimum inhibitory/eradication
concentration (MBIC/MBEC) of unknown antibiofilm peptides. The MBIC model uses an SVM to split the
peptides first, applying an SVR to predicted peptides that are less than 64uM. The MBEC model directly
trains an SVR on the MBEC peptides. Please see the paper for more details.
Both the MBIC and MBEC folders contain source code to train the models which sweep through different
hyperparameters and perform forward selection. Both folders contain the optimized features found during 
forward selection (.json file) and a script that uses the found hyperparameters/features to make predictions 
about unknown peptides. The MBIC folder contains an additional script that combines both the tuned SVM 
and SVR which are then used to perform cross-validation and return an average RMSE of the model. 