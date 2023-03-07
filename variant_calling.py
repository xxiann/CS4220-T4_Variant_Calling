import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split
from sklearn.experimental import enable_halving_search_cv
from sklearn.model_selection import HalvingGridSearchCV
import time

def parse_train(files): # pass sequence of data files
    dataset = []
    for file in files:
        # read VCF data and truth labels
        vcf = pd.read_table('./Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/ZB/CS4220/project_1/data/'+file+'/snv-parse-'+file+'.txt', low_memory=False)
        if file == 'real2_part1':
            truth = pd.read_table('./Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/ZB/CS4220/project_1/data/'+file+'/real2_truth_chr1to5.bed', header=None, usecols=[0, 1], names=['Chr', 'POS'])
        else:
            truth = pd.read_table('./Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/ZB/CS4220/project_1/data/'+file+'/'+file+'_truth.bed', header=None, usecols=[0, 1], names=['Chr', 'POS'])

        # merge VCF file with truth labels
        truth_chr = vcf['Chr'].isin(truth['Chr'])
        truth_pos = vcf['START_POS_REF'].isin(truth['POS'])
        vcf['LABEL'] = pd.concat([truth_chr, truth_pos], axis=1).all(axis=1)

        vcf.drop(columns=['Chr', 'START_POS_REF', 'END_POS_REF', 'REF', 'ALT', 'Sample_Name'], inplace=True)
        dataset.append(vcf)

    train_set = pd.concat(dataset, ignore_index=True)

    # encode categorical variables
    # train_set = pd.get_dummies(train_set, columns=['REF', 'ALT'], prefix=['REF', 'ALT'])

    return train_set

def parse_test(dataset):
    feat = pd.read_table('./Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/ZB/CS4220/project_1/data/'+dataset+'/snv-parse-'+dataset+'.txt', low_memory=False)
    feat_filtered = feat.drop(columns=['Chr', 'START_POS_REF', 'END_POS_REF', 'REF', 'ALT', 'Sample_Name'])
    # feat_encoded = pd.get_dummies(feat_encoded, columns=['REF', 'ALT'], prefix=['REF', 'ALT'])
    return feat, feat_filtered, dataset

def train(feat):
    # choose features to train model on
    X, y = feat.drop('LABEL', axis=1), feat['LABEL'] # use all features except chr no. and genomic coordinate
    # data = xgb.DMatrix(data=X,label=y)
    # X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=123)

    # find optimal parameters using GridSearch
    param_grid = {'n_estimators': [50, 100, 200],
                  "max_depth": [3, 4, 7, 10, 25],
                  "gamma": [0.5, 1, 5, 10, 25],
                  "min_child_weight": [1, 3, 5, 10, 25],
                  "scale_pos_weight": [1, 3, 5, 10, 25]}
    xg_class = xgb.XGBClassifier(gamma=10, max_depth=3, min_child_weight=3, n_estimators=100, scale_pos_weight=5, use_label_encoder=False, eval_metric='logloss') # instantiate classifier
    # xg_class = xgb.XGBClassifier(use_label_encoder=False, eval_metric='logloss')
    # grid_cv = HalvingGridSearchCV(xg_class, param_grid, scoring='roc_auc', n_jobs=-1, verbose=1)
    xg_class.fit(X, y)
    # print(grid_cv.best_params_)
    # print(grid_cv.best_score_)
    # print(xg_class.feature_importances_)
    # print(xg_class.get_booster().get_score(importance_type='gain'))
    xg_class.save_model('./Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/ZB/CS4220/project_1/model.json')

def test(test_set, test_set_filtered, dataset): # input test features
    xg_class = xgb.XGBClassifier(use_label_encoder=False, eval_metric='logloss') # instantiate classifier
    xg_class.load_model('./Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/ZB/CS4220/project_1/model.json')
    labels = xg_class.predict(test_set_filtered) # predicted labels
    # rmse = np.sqrt(mean_squared_error(y_test, preds))
    # print("RMSE: %f" % (rmse))

    preds = pd.DataFrame(columns=['Chr', 'POS'])
    # list out Chr and POS of predicted labels
    for index, label in enumerate(labels):
        if label:
            row = test_set.iloc[index, 0:2]
            row_df = pd.DataFrame({'Chr': row['Chr'], 'POS': row['START_POS_REF']}, index=[0])
            preds = pd.concat([preds, row_df], ignore_index=True)
    preds.to_csv('./Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/ZB/CS4220/project_1/data/'+dataset+'/'+dataset+'_predictions.csv', index=False)

start = time.time()

# train model on a combined train set with syn1-5, real1 and real2_part1
train(parse_train(['syn1', 'syn2', 'syn3', 'syn4', 'syn5', 'real1']))

end = time.time()
elapsed = end - start
print('Time taken for training:', elapsed, 'seconds')

# test model on real2_part2
test(*parse_test('real2_part1'))