#!/usr/bin/env python
# Copyright 2015(c) The Ontario Institute for Cancer Research. All rights reserved.
#
# This program and the accompanying materials are made available under the
# terms of the GNU Public License v3.0.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>.

from __future__ import division
import os, sys
import math
from sklearn import metrics
from sklearn.metrics import *
import numpy as np
import subprocess
from numpy import *
from collections import Counter
from functools import partial
from subprocess import Popen, PIPE, STDOUT



#NORM for ED
def kernel_norm(p,q):
    pq = np.dot(p,q)
    pp = np.dot(p,p.T)
    qq = np.dot(q,q.T)
    norm = pq/(math.sqrt(pp*qq))
    return norm


#ED 
def norm_euclidean(p,q):
    # get frequencies from counts 
    p = np.asarray(p/np.sum(p), dtype=np.float)
    q = np.asarray(q/np.sum(q), dtype=np.float)
    d = kernel_norm(p,p) + kernel_norm(q,q) - 2*kernel_norm(p,q)
    return math.sqrt(d)



# KL sum distance 
def kl_sum_pw(p,q):
    # get rid of zeros in q
    p = p[q != 0]
    q = q[q != 0]
    # get rid of zeros in p
    p = p[p != 0]
    q = q[p != 0]
    #return np.sum(np.where(p != 0, p * np.log(p / q), 0))
    return np.sum(p * np.log(p / q))


# JSD
def js_pairwise(p, q):
    """Jensen-Shannon divergence
    # p and q are k-mer counts 
    Parameters
    ----------
    p, q : array-like, dtype=float, shape=n
    Discrete probability distributions.
    """
    p = np.asarray(p/np.sum(p), dtype=np.float)
    q = np.asarray(q/np.sum(q), dtype=np.float)
    r = (p + q)/2    
    js = 0.5 * kl_sum_pw(p,r) + 0.5 * kl_sum_pw(q,r)
    return js



def get_distance_performance_pairwise(X_train,X_test,y_train,y_test,metric):
    dist = pairwise_distances(X_test, X_train, metric=metric)
    #print dist.shape[0], dist.shape[1], len(y_train), len(np.unique(y_train))
    acc, ftmean = get_accuracy_scores_tune_knn_and_predict(dist,y_train,y_test,K=1) 
    return acc, ftmean



def get_distance_performance_pairwise_loo(X_train,X_test,y_train,y_test,metric):
    dist = pairwise_distances(X_test, X_train, metric=metric)
    #print dist.shape[0], dist.shape[1], len(y_train), len(np.unique(y_train))
    acc, ft = get_accuracy_scores_tune_knn_and_predict_loo(dist,y_train,y_test,K=1) 
    return acc, ft




def get_accuracy_scores_tune_knn_and_predict_loo(distances,y_train,y_test,K):
    n_x = range(0,distances.shape[0])
    f_stat_p = partial(f_stat,Y = y_train)
    for row in n_x:
        v = distances[row,:]
        ft = np.apply_along_axis(f_stat_p,axis=0,arr=array(v))
        y_pred_combined_dist = find_knn(v,y_train,weighted=True,K=K)[0]
    if y_pred_combined_dist == y_test:
        accu = 1
    else:
        accu = 0
    return (accu, ft)



def run_knn_and_predict(distances,y_train, K):
    y_pred_combined_dist=[]
    n_x = range(0,distances.shape[0])
    for row in n_x:
        v = distances[row,:]
        y_pred_combined_dist.append(find_knn(v,y_train,weighted=True,K=K)[0])
    return y_pred_combined_dist



def find_knn(x,classes_train,weighted=False,K=1):
 """ find K nearest neighbours given a distance vector d(train_classes)
 """
 ndata = x.shape[0]
 K = K if K < ndata else ndata
 idx = argsort(x) # sorting
 # return the indexes of K nearest neighbours
 if weighted:
     x[x == 0.0] = 0.000000000000001
     weights = 1.0/x[idx[:K]]
     classes_predicted = array(list(set(classes_train[idx[:K]])))
     # this calculates the weight for each class
     evidence = array([weights[classes_train[idx[:K]] == c].sum() for c in classes_predicted])
     prob_pred = float(evidence.max()/evidence.sum())
     predicted_class = str(classes_predicted[evidence == evidence.max()][0])
 else:
     count = Counter(classes_train[idx[:K]])
     predicted_class = count.most_common(1)[0][0] 
     prob_pred = float(count.most_common(1)[0][1]/K)
 return (predicted_class,prob_pred)



def training_kmer_performance(X_train,y_train,X_test,t_train,t_test):
    distances_train = {}
    dist_accuracies_train = {}
    ftmean_train = {}
    feature = 'counts'
    dist_accuracies_train['d_js'], ftmean_train['d_js'] = get_distance_performance_pairwise(array(X_train[feature])[t_train],array(X_test[feature])[t_test],y_train[t_train],y_train[t_test],metric=js_pairwise)
    dist_accuracies_train['d_eucl'], ftmean_train['d_eucl'] =  get_distance_performance_pairwise(array(X_train[feature])[t_train],array(X_test[feature])[t_test],y_train[t_train],y_train[t_test],metric=norm_euclidean)
    return (dist_accuracies_train, ftmean_train)



def training_kmer_performance_loo(X_train,y_train,X_test,t_train,t_test):
    distances_train = {}
    dist_accuracies_train = {}
    ftmean_train = {}
    feature = 'counts'
    dist_accuracies_train['d_js'], ftmean_train['d_js'] = get_distance_performance_pairwise_loo(array(X_train[feature])[t_train],array(X_test[feature])[t_test],y_train[t_train],y_train[t_test],metric=js_pairwise)
    dist_accuracies_train['d_eucl'], ftmean_train['d_eucl'] =  get_distance_performance_pairwise_loo(array(X_train[feature])[t_train],array(X_test[feature])[t_test],y_train[t_train],y_train[t_test],metric=norm_euclidean)
    return (dist_accuracies_train, ftmean_train)



def return_max_in_array( x ):
    return max(x)



def f_stat(X,Y):
    n = array(X).shape[0]
    classes = np.unique(Y)
    assert(len(Y) == n)
    #alpha=1e-4
    mean_global = np.mean(X, axis=0)
    between_group = 0
    within_group = 0
    for c in classes:
        n_c = len(Y[Y==c])
        if n_c < 2:
            continue
        mu_group = np.mean(X[Y==c], axis=0)
        mu_diff = n_c*(mu_group - mean_global)**2
        between_group = between_group + mu_diff
        within_class = array([(s - mu_group)**2 for s in list(X[Y==c])])
        within_group = within_group + np.sum(within_class)
        if within_group == 0.0:
            f = 0.0
        else:
            f = (between_group/(len(classes) - 1))/(within_group/(n - len(classes))) 
    return f




def knn_predict(X_train, y_train, X_test, dist_blast_test_set, dist_cb_test_set, dist_accuracies_train):
    distances_test = {}
    dist_accuracies_test = {}
    distances_test['d_blast'] = dist_blast_test_set
    distances_test['d_cb'] = dist_cb_test_set
    if X_train != None:
        if 'd_js' in dist_accuracies_train.keys():
            feature = 'counts'
            distances_test['d_js'] = pairwise_distances(array(X_test[feature]), array(X_train[feature]), metric=js_pairwise)
        else:
            distances_test['d_js'] = None
        if 'd_eucl' in dist_accuracies_train.keys():
            feature = 'counts'
            distances_test['d_eucl'] = pairwise_distances(array(X_test[feature]), array(X_train[feature]), metric=norm_euclidean)
        else:
            distances_test['d_eucl'] = None
    kn = 1
    f_stat_p = partial(f_stat,Y = y_train)
    f_sum = 0
    d_combined = 0
    d_count = 0
    for d in distances_test.keys():
        if distances_test[d] is not None:
            d_count += 1   
            f = np.apply_along_axis(f_stat_p,axis=1,arr=distances_test[d])
            f_sum = f_sum + f
            d_combined = d_combined + distances_test[d]*f[:,None]
    if d_count > 0:
        f_sum[f_sum == 0.0] = 1.0
        d_combined = d_combined/f_sum[:,None]
        pred = run_knn_and_predict(d_combined, y_train, K=kn)
    # this in the case all distances have 0 predictiveness
    elif d_count == 0:
        pred = 'NA'
    return (pred, float(f_sum))




def two_sample_t_test_welch(a,b):
    ''' Welch t-test'''
    from scipy.special import stdtr
    abar = a.mean()
    avar = a.var(ddof=1)
    na = a.size
    adof = na - 1
    bbar = b.mean()
    bvar = b.var(ddof=1)
    nb = b.size
    bdof = nb - 1
    # Compute Welch's t-test using the descriptive statistics.
    tf = (abar - bbar) / np.sqrt(avar/na + bvar/nb)
    dof = (avar/na + bvar/nb)**2 / (avar**2/(na**2*adof) + bvar**2/(nb**2*bdof))
    pf = stdtr(dof, -np.abs(tf))
    return (tf, pf)



def compress_data_train_mp(seq2,dir_out,plzip_threads):
        '''
        Compress the training set 
        '''
        seq = str(seq2.seq)
        gi_id2 = str(seq2.id.split("|")[1])
        p = Popen(['plzip', '-n', str(plzip_threads), '-4'], stdout=PIPE, stdin=PIPE, stderr=STDOUT)    
        len_compressed = len(p.communicate(input=seq)[0])
        return (gi_id2, len_compressed) 



def get_ncd_distance_mp_train(seq2, seq1, plzip_threads, dict_train):
    '''
    Get the compression scores
    '''
    gi_id2 = str(seq2.id.split("|")[1])
    gi_id1 = str(seq1.id.split("|")[1])
    cy = dict_train[gi_id2]
    cx = dict_train[gi_id1]
    sc = str(seq1.seq) + str(seq2.seq) 
    p = Popen(['plzip', '-n', str(plzip_threads), '-4'], stdout=PIPE, stdin=PIPE, stderr=STDOUT)    
    cxy = len(p.communicate(input=sc)[0])
    if cy > cx:
        n = (cxy - cx) / cy
    else:
        n = (cxy - cy) / cx
    return n
    


def get_ncd_distance_mp_test(seq2, seq1, plzip_threads, dict_train):
    gi_id_train = str(seq2.id.split("|")[1])
    cy = dict_train[gi_id_train]
    sc = str(seq1.seq) + str(seq2.seq)
    p = Popen(['plzip', '-n', str(plzip_threads), '-4'], stdout=PIPE, stdin=PIPE, stderr=STDOUT)    
    cxy = len(p.communicate(input=sc)[0])
    p = Popen(['plzip', '-n', str(plzip_threads), '-4'], stdout=PIPE, stdin=PIPE, stderr=STDOUT)    
    s_test = str(seq1.seq)
    cx = len(p.communicate(input=s_test)[0])
    if cy > cx:
       n = (cxy - cx) / cy
    else:
       n = (cxy - cy) / cx
    return n


