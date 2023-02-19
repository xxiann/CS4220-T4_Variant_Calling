import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split

def parse_train(dataset):
    # read VCF data and truth labels
    vcf = pd.read_table('./Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/ZB/CS4220/project_1/data/'+dataset+'/snv-parse-'+dataset+'.txt', low_memory=False)
    truth = pd.read_table('./Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/ZB/CS4220/project_1/data/'+dataset+'/'+dataset+'_truth.bed', header=None, usecols=[0, 1], names=['Chr', 'POS'])
    vcf.drop(columns=['END_POS_REF', 'Sample_Name'], inplace=True)
    vcf.rename(columns={'START_POS_REF': 'POS'}, inplace=True)

    # merge VCF file with truth labels
    truth_chr = vcf['Chr'].isin(truth['Chr'])
    truth_pos = vcf['POS'].isin(truth['POS'])
    vcf['LABEL'] = pd.concat([truth_chr, truth_pos], axis=1).all(axis=1)

    # encode categorical variables
    vcf = pd.get_dummies(vcf, columns=['Chr', 'REF', 'ALT'], prefix=['Chr', 'REF', 'ALT'])
    return vcf

def parse_test(dataset):
    feat = pd.read_table('./Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/ZB/CS4220/project_1/data/'+dataset+'/snv-parse-'+dataset+'.txt', low_memory=False)
    feat.drop(columns=['END_POS_REF', 'Sample_Name'], inplace=True)
    feat.rename(columns={'START_POS_REF': 'POS'}, inplace=True)
    feat_encoded = pd.get_dummies(feat, columns=['Chr', 'REF', 'ALT'], prefix=['Chr', 'REF', 'ALT'])
    return feat, feat_encoded, dataset

def train(feat, existing_model):
    # choose features to train model on
    X, y = feat.drop('LABEL', axis=1), feat['LABEL'] # use all features
    # data = xgb.DMatrix(data=X,label=y)
    # X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=123)
    xg_class = xgb.XGBClassifier(use_label_encoder=False) # instantiate classifier
    if existing_model:
        xg_class.fit(X, y, xgb_model='./Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/ZB/CS4220/project_1/model.json')
    else:
         xg_class.fit(X, y)
    xg_class.save_model('./Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/ZB/CS4220/project_1/model.json')

def test(test_set, test_set_encoded, dataset): # input test features
    xg_class = xgb.XGBClassifier(use_label_encoder=False) # instantiate classifier
    xg_class.load_model('./Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/ZB/CS4220/project_1/model.json')
    # print(test_set.shape)
    labels = xg_class.predict(test_set_encoded) # predicted labels
    # rmse = np.sqrt(mean_squared_error(y_test, preds))
    # print("RMSE: %f" % (rmse))

    preds = pd.DataFrame(columns=['Chr', 'POS'])
    # list out Chr and POS of predicted labels
    for index, label in enumerate(labels):
        if label:
            row = test_set.iloc[index, 0:2]
            row_df = pd.DataFrame({'Chr': row['Chr'], 'POS': row['POS']}, index=[0])
            preds = pd.concat([preds, row_df], ignore_index=True)
    # print(preds)
    preds.to_csv('./Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/ZB/CS4220/project_1/data/'+dataset+'/'+dataset+'_predictions.csv', index=False)

# train model on syn2
train(parse_train('syn2'), False)
print("trained on syn2")

# train model on syn2 to syn4
# for i in range(2, 5):
#     train(parse_train('syn'+str(i)), True)
#     print("trained on syn"+str(i))

# test model on syn3
test(*parse_test('syn3'))

# deal with features which are present in some training data only / present in test data only