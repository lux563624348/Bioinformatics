import numpy as np
import pandas as pd
from pandas.tseries.offsets import *
import os
import re
import argparse
import time
import itertools
import sys

import sklearn
from sklearn.preprocessing import LabelEncoder

from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.feature_extraction import DictVectorizer
from sklearn.preprocessing import StandardScaler
from sklearn.feature_extraction import FeatureHasher
from sklearn.feature_selection import SelectKBest

from sklearn.linear_model import LogisticRegression

from scipy.special import logit
from scipy.special import expit

from minepy import MINE

import xgboost as xgb

import scipy.sparse
import scipy.io


from imblearn.under_sampling import RandomUnderSampler
from imblearn.ensemble import EasyEnsemble
from imblearn.ensemble import BalanceCascade
    
def label_encoding(df, cols):
    for col in cols:
        le = LabelEncoder()
        le.fit(df[col])
        col_numerical_data = le.transform(df[col])
        df = df.drop(col, axis=1)
        df[col] = col_numerical_data
    return df

def onehot_encoding(df, cols):
    vec = DictVectorizer()
    
    vec_data = pd.DataFrame(vec.fit_transform(df[cols].astype(str).to_dict(orient='records')).toarray())
    vec_data.columns = vec.get_feature_names()
    vec_data.index = df.index
    
    df = df.drop(cols, axis=1)
    df = df.join(vec_data)
    return df 

def create_match_feature(df, col1, col2, new_col):
    df[new_col] = (df[col1] == df[col2]).astype(int)
    return None

def create_diff_feature(df, col1, col2, new_col):
    df[new_col] = df[col1] - df[col2]
    return None

def create_poly2_feature(df, col1, col2, new_col='default'):
    if new_col == 'default':
        new_col = col1 + ',' + col2
    df[new_col] = df.apply(lambda row: str(row[col1]) + ',' + str(row[col2]), axis=1)
    return None

def standarize_feature(df, cols):
    scaler = StandardScaler()
    for col in cols:
        df[col] = scaler.fit_transform(df[col].reshape(-1,1).astype(np.float32))
    return None

def extend_bounds(bins):
    bins[0] = bins[0] - 1
    bins[-1] = bins[-1] + 1

def discretize_feature(df, cols, num_bins, how='equal_freq'):
    if how == 'equal_width':
        for col in cols:
            group_names = range(num_bins)
            df[col] = pd.cut(df[col], num_bins, labels=group_names)
    elif how == 'equal_freq':
        for col in cols:
            freq_per_bin = df[col].shape[0] / num_bins
            values = sorted(df[col])
            bins = [values[0]]
            group_names = [0]
            for i in range(1, num_bins + 1):
                if i < num_bins:
                    bins.append(values[i * freq_per_bin - 1])
                else:
                    bins.append(values[-1])
            group_names = range(num_bins)
         
            i = 1
            while i < len(bins):
                if bins[i] == bins[i-1]:
                    del bins[i]
                    del group_names[i]
                else:
                    i += 1
            extend_bounds(bins)
            df[col] = pd.cut(df[col], bins, labels=group_names)
    return None


def load_data(subsample=1, random=True, length=False, distance_to_TES=False, kmer_freq=False, splice_site_score=False, RBP_motif=False, RBP_binding_score_eclip=False, polII=False):
    data_path = '/home/lxiang/cloud_research/PengGroup/ZZeng/Data/intron_retention/data/Human_summary/IRI/introns'
    
    num_samples = len(os.listdir(data_path))
    df_list = []
    for data_full_path in [os.path.join(data_path, f) for f in os.listdir(data_path) if re.search('introns.txt$', f)]:
        df = pd.read_csv(data_full_path, header=0, sep='\t', na_values=['NA', "NA (5'AS)", "NA (3'AS)", "NA (unannotated exon)"]).loc[:,['CIR_id', 'adjacent_CER_RPKM', 'intron_IRI']]
        df = df[(df.intron_IRI >= 0) & (df.intron_IRI <= 1) & (df.adjacent_CER_RPKM >= 0.1)].drop(['adjacent_CER_RPKM'], axis=1).set_index('CIR_id')
        df_list.append(df)    
        
    IRI_df = pd.concat(df_list, axis=1).dropna()
    IRI_df['mean_intron_IRI'] = IRI_df.apply(lambda row: np.mean(row[:num_samples]), axis=1)
    
    
    num_IRI_high_CIRs = 5000
    num_IRI_low_CIRs = 5000
    
    IRI_high_CIR_list = IRI_df.sort_values(by='mean_intron_IRI', ascending=False).iloc[:num_IRI_high_CIRs].index.tolist()
    IRI_low_CIR_list = IRI_df.sort_values(by='mean_intron_IRI').iloc[:num_IRI_low_CIRs].index.tolist()
    
    train_test_df = pd.DataFrame({'CIR_id': IRI_high_CIR_list + IRI_low_CIR_list, 'intron_IRI': [1] * num_IRI_high_CIRs + [0] * num_IRI_low_CIRs})
    train_test_df['gene_id'] = train_test_df.CIR_id.str[:-4]
    
    # length
    if length:
        # gene length
        genes_length_path = '/home/lxiang/cloud_research/PengGroup/ZZeng/Data/intron_retention/ML_data/data/IRstat.length.genes.txt'
        genes_length_df = pd.read_csv(genes_length_path, sep='\t', header=0)
    
        # CIR length
        introns_length_path = '/home/lxiang/cloud_research/PengGroup/ZZeng/Data/intron_retention/ML_data/data/IRstat.length.introns.txt'
        introns_length_df = pd.read_csv(introns_length_path, sep='\t', header=0)  
        
        train_test_df = train_test_df.merge(genes_length_df, on='gene_id', how='left')
        train_test_df = train_test_df.merge(introns_length_df, on='CIR_id', how='left')                
    
    # distance to TES
    if distance_to_TES:
        distance_to_TES_path = '/home/lxiang/cloud_research/PengGroup/ZZeng/Data/intron_retention/ML_data/data/IRstat.distance_to_TES.introns.txt'
        distance_to_TES_df = pd.read_csv(distance_to_TES_path, sep='\t', header=0)
        
        train_test_df = train_test_df.merge(distance_to_TES_df, on='CIR_id', how='left')
        
    # kmer freq (1-mers, 2-mers, 3-mers)
    if kmer_freq:
        introns_kmer_freq_path = '/home/lxiang/cloud_research/PengGroup/ZZeng/Data/intron_retention/ML_data/data/IRstat.kmer_freq.introns.txt'
        introns_kmer_freq_df = pd.read_csv(introns_kmer_freq_path, sep='\t', header=0)   
        
        train_test_df = train_test_df.merge(introns_kmer_freq_df, on='CIR_id', how='left')
        
    # splice site score
    if splice_site_score:
        introns_splice_site_score_path = '/home/lxiang/cloud_research/PengGroup/ZZeng/Data/intron_retention/ML_data/data/IRstat.splice_site_score.introns.txt'
        introns_splice_site_score_df = pd.read_csv(introns_splice_site_score_path, sep='\t', header=0)
        
        train_test_df = train_test_df.merge(introns_splice_site_score_df.loc[:,['CIR_id', '5_splice_site_score', '3_splice_site_score']], on='CIR_id', how='left')
    
    # RBP motif    
    if RBP_motif:
        introns_RBP_motif_upstream_CIR_path = '/home/lxiang/cloud_research/PengGroup/ZZeng/Data/intron_retention/data/exploration/motif/fimo/IRstat.RBP_motif.upstream_CIR.introns.txt'
        introns_RBP_motif_upstream_CIR_df = pd.read_csv(introns_RBP_motif_upstream_CIR_path, sep='\t', header=0)
        introns_RBP_motif_upstream_CIR_df.rename(columns=dict(zip(introns_RBP_motif_upstream_CIR_df.columns[1:], ['upstream_CIR_'+x for x in introns_RBP_motif_upstream_CIR_df.columns[1:]])), inplace=True)
        
        introns_RBP_motif_downstream_CIR_path = '/home/lxiang/cloud_research/PengGroup/ZZeng/Data/intron_retention/data/exploration/motif/fimo/IRstat.RBP_motif.downstream_CIR.introns.txt'
        introns_RBP_motif_downstream_CIR_df = pd.read_csv(introns_RBP_motif_downstream_CIR_path, sep='\t', header=0)
        introns_RBP_motif_downstream_CIR_df.rename(columns=dict(zip(introns_RBP_motif_downstream_CIR_df.columns[1:], ['downstream_CIR_'+x for x in introns_RBP_motif_downstream_CIR_df.columns[1:]])), inplace=True)    
        
        train_test_df = train_test_df.merge(introns_RBP_motif_upstream_CIR_df, on='CIR_id', how='left')
        train_test_df = train_test_df.merge(introns_RBP_motif_downstream_CIR_df, on='CIR_id', how='left')
    
    # RBP motif binding score eclip   
    if RBP_binding_score_eclip:      
        introns_RBP_binding_score_eclip_path = '/home/lxiang/cloud_research/PengGroup/ZZeng/Data/intron_retention/ML_data/data/IRstat.RBP_binding_score_eclip_K562.introns.txt'
        introns_RBP_binding_score_eclip_df = pd.read_csv(introns_RBP_binding_score_eclip_path, sep='\t', header=0)   
        
        train_test_df = train_test_df.merge(introns_RBP_binding_score_eclip_df, on='CIR_id', how='left')  
        
    # polII
    if polII:
        introns_polII_path = '/home/lxiang/cloud_research/PengGroup/ZZeng/Data/intron_retention/ML_data/data/Rest_CD4_polll.quant.IRI.introns.txt'
        introns_polII_df = pd.read_csv(introns_polII_path, header=0, sep='\t').loc[:,['CIR_id', 'CIR_RPKM', 'intron_IRI']].rename(columns={'CIR_RPKM': 'pollII_RPKM', 'intron_IRI': 'polII_IRI'})
        
        train_test_df = train_test_df.merge(introns_polII_df, on='CIR_id', how='left') 
     
    train_test_df.drop(labels='gene_id', axis=1, inplace=True)    
        
    if random:         
        subsample_indice = np.random.choice(train_test_df.shape[0], int(train_test_df.shape[0] * subsample), replace=False)
    else:
        subsample_indice = range(int(train_test_df.shape[0] * subsample))
    train_test_df = train_test_df.iloc[subsample_indice,:]
    
    return train_test_df
           
    
def process_train_test_df(train_test_df, train_size=0.8, encoding='onehot_encoding', remove_missing_rows=True, fill_missing_values=False, discretization=False, standarization=False, poly2_feature=False, select_features=False):
    num_train = int(train_test_df.shape[0] * train_size)
    
    if remove_missing_rows:
        num_train = train_test_df.iloc[:num_train,:].replace({np.inf: np.nan}).dropna().shape[0]
        train_test_df = train_test_df.replace({np.inf: np.nan}).dropna()
        
    if fill_missing_values:
        train_test_df = train_test_df.apply(lambda col: col.fillna(col.value_counts().index[0]), axis=0)
    
    train_indice = np.array(train_test_df['CIR_id'][:num_train])
    test_indice = np.array(train_test_df['CIR_id'][num_train:])
    
    y = np.array(train_test_df['intron_IRI'])
    train_test_df = train_test_df.drop(['intron_IRI', 'CIR_id'], axis=1)
        
    numerical_features = train_test_df.columns
    categorical_features = [item for item in train_test_df.columns if item not in numerical_features]
    
    if standarization:
        standarized_features = numerical_features
        standarize_feature(train_test_df, standarized_features)    
            
    if discretization:
        discretized_features = numerical_features
        discretize_feature(train_test_df, discretized_features, num_bins=10, how='equal_freq')  
        
    if poly2_feature:
        comb_feature_list = train_test_df.columns
        comb_feature_pairs = list(itertools.combinations(comb_feature_list, 2))
        random_indice = np.random.choice(range(len(comb_feature_pairs)), 20, replace=False)
        random_comb_feature_pairs = [comb_feature_pairs[i] for i in random_indice]
     
        for col1, col2 in random_comb_feature_pairs:
            create_poly2_feature(train_test_df, col1, col2)
            categorical_features.append(col1 + ',' + col2)
       
    if encoding == 'label_encoding':
        train_test_df = label_encoding(train_test_df, categorical_features)      
    
    elif encoding == 'onehot_encoding':
        train_test_df = onehot_encoding(train_test_df, categorical_features)
    
    X = np.array(train_test_df)
    
    if select_features:
        X, feature_list = feature_selection(X, y, num_train, list(train_test_df.columns))
        
    X_train = X[:num_train, :]
    y_train = y[:num_train]

    X_test = X[num_train:, :]
    y_test = y[num_train:]
    
    if select_features:
        return X_train, X_test, y_train, y_test, feature_list, train_indice, test_indice
    else:
        return X_train, X_test, y_train, y_test, list(train_test_df.columns), train_indice, test_indice       

def mic(x, y):
    m = MINE()
    m.compute_score(x, y)
    return (m.mic(), 0.5)

def xgboost_eval(X, y, num_train, indice_list):
    dtrain = xgb.DMatrix(X[:num_train,indice_list], label=y[:num_train], missing=np.nan)
    dtest = xgb.DMatrix(X[num_train:,indice_list], label=y[num_train:], missing=np.nan)
    num_rounds = 30
    progress = {}
    watchlist = [(dtest,'test')]
    params = {"objective": 'binary:logistic',
              "booster" : "gbtree",        
              "silent": 1,
              "seed": 0,
              "nthread": 8,
              "eval_metric": ['auc']
              } 
    model = xgb.train(params, dtrain, num_boost_round=num_rounds, verbose_eval=False, evals=watchlist, evals_result=progress)
    return float(progress['test']['auc'][-1])
    

def feature_selection(X, y, num_train, feature_list, MI_filter = 100, forward_selection = 20, backward_filter = None, xgb_feature_importance_filter = None, markov_blanket_filter = None):
    print "num of features before selection: %d" % len(feature_list)
    if MI_filter:
        print "start MI (mutual information) filter:"
        MI_fit = SelectKBest(lambda X, Y: np.array(map(lambda x:mic(x, Y), X.T)).T, k=len(feature_list) - MI_filter)
        MI_fit.fit(X[:num_train],y[:num_train])
        remain_indice_list = [indice for indice in range(len(feature_list)) if indice in MI_fit.scores_.argsort()[MI_filter:]]
        filter_feature_list = [feature_list[indice] for indice in range(len(feature_list)) if indice not in remain_indice_list]
        feature_list = [feature_list[indice] for indice in remain_indice_list]
        X = X[:,remain_indice_list]
        print "%d features filtered: %s" % (MI_filter, ' '.join(filter_feature_list))
        sys.stdout.flush()
        
        if not forward_selection and not backward_filter:
            return X, feature_list
    
    if forward_selection:
        print "start forward selection:"
        sys.stdout.flush()
        early_stop = False
        select_indice_list = []
        remain_indice_list = range(len(feature_list))
        eval_list = []
        
        for i in range(forward_selection):
            best_indice = None
            if not early_stop or i == 0:
                best_eval = 0
            else:
                best_eval = eval_list[-1] 
                
            for test_indice in remain_indice_list:
                test_eval = xgboost_eval(X, y, num_train, select_indice_list + [test_indice])
                if test_eval > best_eval:
                    best_indice = test_indice
                    best_eval = test_eval
                sys.stdout.flush()
            
            if early_stop and i != 0 and eval_list[-1] >= best_eval:
                print "This iteration doesn't improve test AUC. Early stop!"
                print "%d features selected: %s" % (len(select_indice_list), ' '.join([feature_list[indice] for indice in select_indice_list]))
                sys.stdout.flush()  
                return X[:,select_indice_list], [feature_list[indice] for indice in select_indice_list]
            
            print "%d / %d" % (i+1, forward_selection)
            print "best feature: %s" % feature_list[best_indice]
            print "best test AUC: %f" % best_eval
            sys.stdout.flush()            
            
            eval_list.append(best_eval)
            select_indice_list.append(best_indice)
            remain_indice_list.remove(best_indice)
    
        print "%d features selected: %s" % (len(select_indice_list), ' '.join([feature_list[indice] for indice in select_indice_list]))
        sys.stdout.flush() 
        return X[:,select_indice_list], [feature_list[indice] for indice in select_indice_list]
    
    if backward_filter:
        print "start backward filter:"
        sys.stdout.flush()
        early_stop = True        
        remain_indice_list = range(len(feature_list))
        eval_list = []
      
        for i in range(backward_filter):
            worst_indice = None
            if not early_stop or i == 0:
                best_eval = 0
            else:
                best_eval = eval_list[-1]      
                
            for test_indice in remain_indice_list:
                test_eval = xgboost_eval(X, y, num_train, [indice for indice in remain_indice_list if indice != test_indice])
                if test_eval > best_eval:
                    worst_indice = test_indice
                    best_eval = test_eval
            
            if early_stop and i != 0 and eval_list[-1] >= best_eval:
                print "This iteration doesn't improve test AUC. Early stop!"
                print "%d features remained: %s" % (len(remain_indice_list), ' '.join([feature_list[indice] for indice in remain_indice_list]))
                sys.stdout.flush()  
                return X[:,remain_indice_list], [feature_list[indice] for indice in remain_indice_list]  
            
            print "%d / %d" % (i+1, backward_filter)
            print "worst feature: %s" % feature_list[worst_indice]
            print "best test AUC if remove worst feature: %f" % best_eval
            sys.stdout.flush()            

            eval_list.append(best_eval)
            remain_indice_list.remove(worst_indice)
    
        print "%d features remained: %s" % (len(remain_indice_list), ' '.join([feature_list[indice] for indice in remain_indice_list]))   
        sys.stdout.flush() 
        return X[:,remain_indice_list], [feature_list[indice] for indice in remain_indice_list]
    
    if xgb_feature_importance_filter:
        print "start xgb feature importance filter:"
        sys.stdout.flush()        
        xgb_feature_importance_list =  []
        filter_feature_list = xgb_feature_importance_list[-xgb_feature_importance_filter:]
        filter_indice_list = [feature_list.index(feat) for feat in filter_feature_list]
        remain_indice_list = [indice for indice in range(len(feature_list)) if indice not in filter_indice_list]
        
        print "%d features filtered: %s" % (xgb_feature_importance_filter, ' '.join(filter_feature_list))
        sys.stdout.flush() 
        return X[:,remain_indice_list], [feature_list[indice] for indice in remain_indice_list]
    
    if markov_blanket_filter:
        print "start markov blanket filter:"
        sys.stdout.flush()          
        filter_feature_list = []
        filter_indice_list = [feature_list.index(feat) for feat in filter_feature_list]
        remain_indice_list = [indice for indice in range(len(feature_list)) if indice not in filter_indice_list]
        
        print "%d features filtered: %s" % (len(filter_feature_list), ' '.join(filter_feature_list))
        sys.stdout.flush() 
        return X[:,remain_indice_list], [feature_list[indice] for indice in remain_indice_list]        
          
def create_feature_map(feature_list, fmap_name):
    with open(fmap_name, 'w') as f:
        for i, feat in enumerate(feature_list):
            f.write('{0}\t{1}\tq\n'.format(i, feat))

def elapsed_time(start_time, end_time):
    elapsed_sec = end_time - start_time
    h = int(elapsed_sec / (60 * 60))
    m = int((elapsed_sec % (60 * 60)) / 60)
    s = int(elapsed_sec % 60)
    return "{}:{:>02}:{:>02}".format(h, m, s)

def train(X_train, y_train, X_test, y_test, model_type='xgboost:gblinear', verbose_eval=False):
    if model_type == 'LR':
        model = LogisticRegression(n_jobs=-1, penalty='l1')
        model.fit(X_train, y_train)
        
        y_train_pred_prob = model.predict_proba(X_train)[:,1] 
        print "train auc:", roc_auc_score(y_train, y_train_pred_prob)
        print
        #print classification_report(y_train, model.predict_proba(X_train)[:,1]  >= 0.5)    
    else:
        if model_type == 'xgboost:gblinear':
            params = {"objective": 'binary:logistic',
                      "booster" : "gblinear",        
                      "silent": 1,
                      "seed": 0,
                      "nthread": 8,
                      "eval_metric": ['logloss', 'auc'],
                      "alpha": 10
                      }
        elif model_type == 'xgboost:gbtree':
            params = {"objective": 'binary:logistic',
                      "booster" : "gbtree",        
                      "silent": 1,
                      "seed": 0,
                      "nthread": 8,
                      "eval_metric": ['logloss', 'auc']
                      }
            
        dtrain = xgb.DMatrix(X_train, label=y_train, missing=np.nan)    
        dtest = xgb.DMatrix(X_test, label=y_test, missing=np.nan)        
      
        num_rounds = 1000
        early_stopping_rounds = 20
        progress = {}
        watchlist = [(dtrain,'train'), (dtest,'test')]
    
        model = xgb.train(params, dtrain, num_boost_round=num_rounds, verbose_eval=verbose_eval, evals=watchlist, evals_result=progress, early_stopping_rounds=early_stopping_rounds)
        
        y_train_pred_prob = model.predict(dtrain) 
        print "train auc:", roc_auc_score(y_train, y_train_pred_prob)
        print          
        
    return model

def predict(model, model_type, X_test, y_test, output_auc=True):
    if model_type == 'LR':
        y_test_pred_prob = model.predict_proba(X_test)[:,1]      
    else:
        dtest = xgb.DMatrix(X_test, label=y_test, missing=np.nan)   
        y_test_pred_prob = model.predict(dtest)  
        
    if output_auc:
        print "test auc:", roc_auc_score(y_test, y_test_pred_prob)
        print
        
    return y_test_pred_prob
    
def random_downsample(X_train, y_train, pos_neg_ratio=1):
    sampler = RandomUnderSampler(ratio=pos_neg_ratio, random_state=0)
    X_train, y_train = sampler.fit_sample(X_train, y_train)   
    return X_train, y_train

def easy_ensemble(X_train, y_train, X_test, y_test, model_type='xgboost:gblinear', n_subsets=None):
    if n_subsets is None:
        num_pos = np.sum(y_train == 1)
        num_neg = np.sum(y_train == 0)
        n_subsets = num_neg / num_pos
    sampler = EasyEnsemble(n_subsets=n_subsets)
    X_train_list, y_train_list = sampler.fit_sample(X_train, y_train)
    model_list = []
    for i in range(n_subsets):
        model = train(X_train_list[i], y_train_list[i], X_test, y_test, model_type=model_type)
        model_list.append(model)
        
    return model_list
        
def balance_cascade(X_train, y_train, X_test, y_test, model_type='xgboost:gblinear', n_max_subset=None, classifier='gradient-boosting'):
    sampler = BalanceCascade(n_max_subset=n_max_subset, classifier=classifier)
    X_train_list, y_train_list = sampler.fit_sample(X_train, y_train)
    n_subsets = len(y_train_list)
    model_list = []
    for i in range(n_subsets):
        model = train(X_train_list[i], y_train_list[i], X_test, y_test, model_type=model_type)
        model_list.append(model)
        
    return model_list

def ensemble_predict(model_list, model_type, X_test, y_test, output_auc=True):
    if model_type == 'LR':
        y_test_pred_prob = np.mean(map(lambda model: model.predict_proba(X_test)[:,1], model_list), axis=0)
        #y_test_pred_prob = expit(np.mean(map(lambda model: logit(model.predict_proba(X_test)[:,1]), model_list), axis=0))
    else:
        dtest = xgb.DMatrix(X_test, label=y_test, missing=np.nan)
        y_test_pred_prob = np.mean(map(lambda model: model.predict(dtest), model_list), axis=0)
        #y_test_pred_prob = expit(np.mean(map(lambda model: logit(model.predict(dtest)), model_list), axis=0))
        
    if output_auc:               
        print "test auc:", roc_auc_score(y_test, y_test_pred_prob)
        print
    
    return y_test_pred_prob

def find_min_loss_threshold(y_val, yval_pred_prob, fp_loss=1, fn_loss=10):
    fpr, tpr, thresholds = roc_curve(y_val, yval_pred_prob, pos_label=1)
    num_pos = np.sum(y_val == 1)
    num_neg = np.sum(y_val == 0)
    threshold = thresholds[np.argmin((1 - tpr) * num_pos * fn_loss + fpr * num_neg * fp_loss)]
    return threshold

def output_feature_importance(model, model_type, feature_list):
    if model_type == 'LR':
        weight = zip(feature_list, model.coef_[0])
        sorted_weight = sorted(weight, key=lambda x: x[1], reverse=True)
        print "features with positive weight:"
        print [item for item in sorted_weight if item[1] > 0]
        print
        print "features with negative weight:"
        print [item for item in sorted_weight if item[1] < 0][::-1]
        print
        print "features with 0 weight (no contribution):"
        print [item for item in sorted_weight if item[1] == 0]               
        print
        
    elif model_type == 'xgboost:gblinear':
        model_name = 'xgboost.txt'
        model.dump_model(model_name)
        with open(model_name, 'r') as f:
            weight_list = [float(w) for w in re.search('weight:(.+)', f.read(), flags=re.DOTALL).group(1).strip().split('\n')]
        weight = zip(feature_list, weight_list)
        sorted_weight = sorted(weight, key=lambda x: x[1], reverse=True)
        os.remove(model_name)
        print "features with positive weight:"
        print [item for item in sorted_weight if item[1] > 0]
        print
        print "features with negative weight:"
        print [item for item in sorted_weight if item[1] < 0][::-1]
        print
        print "features with 0 weight (no contribution):"
        print [item for item in sorted_weight if item[1] == 0]               
        print   
        
    elif model_type == 'xgboost:gbtree':
        fmap_name = 'xgboost.map'
        create_feature_map(feature_list, fmap_name)
        importance = model.get_fscore(fmap=fmap_name)
        norm_importance = [(item[0], item[1] * 1.0 / sum(importance.values())) for item in importance.items()]
        sorted_norm_importance = sorted(norm_importance, key=lambda x: x[1], reverse=True)
        os.remove(fmap_name)
        print "feature importances:"
        print sorted_norm_importance
        print
           
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--model_name', type=str, help='model saved name')
    parser.add_argument('--model_type', type=str, help='xgboost:gblinear, xgboost:gbtree, LR', default='xgboost:gblinear')
    parser.add_argument('--downsample_type', type=str, help='random_downsample:ratio, easy_ensemble, balance_cascade', default='random_downsample:1')
    parser.add_argument('--threshold', action="store_true", help='determine classification threshold')
    parser.add_argument('--prediction_result', type=str, help='output prediction result')
    parser.add_argument('--seed', type=int, help='random seed', default=0)
    parser.add_argument('--features', type=str, help='i.e., kmer_freq,RBP_binding_score_eclip')
    args = parser.parse_args()   
    
    np.random.seed(args.seed)
    
    start_time = time.time()
    
    model_type = args.model_type
    
    if args.features:
        feature_list = args.features.split(',')
        feature_dict = dict([(feature, True) for feature in feature_list])
        train_test_df = load_data(**feature_dict)
    else:
        train_test_df = load_data()
    
    if args.model_type == 'LR':
        X_train, X_test, y_train, y_test, feature_list, train_indice, test_indice = process_train_test_df(train_test_df, standarization=True)
    else: 
        X_train, X_test, y_train, y_test, feature_list, train_indice, test_indice = process_train_test_df(train_test_df)
        
    print X_train.shape, X_test.shape
    
    if args.threshold:
        X_train,  X_val, y_train, y_val = train_test_split(X_train, y_train, test_size=0.1)
    
    if args.downsample_type:
        downsample_type = args.downsample_type
        if re.match('random_downsample', downsample_type):
            neg_pos_ratio = float(re.search(':(\S+)$', downsample_type).group(1))
            pos_neg_ratio = 1.0 / neg_pos_ratio
            X_train, y_train = random_downsample(X_train, y_train, pos_neg_ratio=pos_neg_ratio)
            model = train(X_train, y_train, X_test, y_test, model_type=model_type)
            y_test_pred_prob = predict(model, model_type, X_test, y_test, output_auc=True)
            output_feature_importance(model, model_type, feature_list)
        elif downsample_type == 'easy_ensemble':
            model_list = easy_ensemble(X_train, y_train, X_test, y_test, model_type=model_type, n_subsets = None)
            y_test_pred_prob = ensemble_predict(model_list, model_type, X_test, y_test)
        elif downsample_type == 'balance_cascade':
            model_list = balance_cascade(X_train, y_train, X_test, y_test, model_type=model_type, n_max_subset=10)
            y_test_pred = ensemble_predict(model_list, model_type, X_test, y_test)
        else:
            print "Wrong value for downsample_type. Exit!"
            sys.exit(1)
    else:
        model = train(X_train, y_train, X_test, y_test, model_type=model_type)
        y_test_pred_prob = predict(model, model_type, X_test, y_test, output_auc=True)
        output_feature_importance(model, model_type, feature_list)        
    
    if args.threshold:
        fp_loss = 1
        fn_loss = 5
        if args.random_downsample:
            yval_pred_prob = predict(model, model_type, X_val, y_val, output_auc=False)
        elif args.easy_ensemble or args.balance_cascade:
            yval_pred_prob = ensemble_predict(model_list, model_type, X_val, y_val, output_auc=False)
        
        threshold = find_min_loss_threshold(y_val, yval_pred_prob, fp_loss=fp_loss, fn_loss=fn_loss)
        
        print "optimal threshold:", threshold
        print 
    
        y_test_pred = (y_test_pred_prob >= threshold).astype(int) 
            
        print "test evaluation:"
        print classification_report(y_test, y_test_pred)
        print
        
    if args.model_name:
        model.save_model(args.model_name)
        print "model params:"
        print params
        print
        print
    
    if args.prediction_result:
        prediction_df = pd.DataFrame({'indice': test_indice, 'pred': y_test_pred_prob, 'true_label': y_test})
        prediction_df.to_csv(args.prediction_result, index=None, header=None, sep='\t')
    
    end_time = time.time()
    #print "Run complete: %s elapsed" % elapsed_time(start_time, end_time)
    print
    print
    
if __name__ == '__main__':
    main()
    
    
