#!/usr/bin/env python3

#############
# Libraries #
#############

import os
import datetime
import argparse
#Suppress lots of memory warning
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 
from os.path import expanduser as eu
from os import path
from math import ceil
import tensorflow as tf
import pandas as pd
import numpy as np
from scipy.stats import poisson 
from scipy.special import loggamma
from scipy.io import mmread
import scipy.sparse as sparse
import warnings
tf.compat.v1.disable_eager_execution() 

from sklearn.decomposition import PCA

def run_pca(df_ref, npcs=1, class_key='gastr_type', class_values=None, timepoint_key='timepoint', timepoint_values=None, celltype_key='cluster', individual_key='MULTI_class') :
    # subset to timepoint and gastruloid type and calculate cell numbers
    tp_ref = df_ref
    if class_values is not None:
        tp_ref = tp_ref[tp_ref[class_key].isin(class_values)]
    if class_values is not None:
        tp_ref = tp_ref[tp_ref[timepoint_key].isin(timepoint_values)]
    tp_ref = tp_ref.value_counts([celltype_key, individual_key]).reset_index(name='counts').pivot_table(index=individual_key, columns=celltype_key, values='counts', fill_value=0)
    for clu in np.unique(df_ref[celltype_key]):
        if clu not in tp_ref:
            tp_ref[clu] = 0
    tp_ref = tp_ref[np.unique(df_ref[celltype_key])]
    tp_ref_vals = np.array(tp_ref)/np.sum(np.array(tp_ref), axis=1).reshape(-1,1)
    tp_ref_vals[tp_ref_vals==0] = min(tp_ref_vals[tp_ref_vals!=0])#*0.1
    tp_ref_vals = np.log(tp_ref_vals)
    
    # prep for PCA
    offset = np.mean(tp_ref_vals, axis=0)
    tp_ref_vals = (tp_ref_vals.T - offset.reshape(-1,1)).T
    
    # perform PCA
    if int(npcs) > 0:
        pca = PCA(n_components=int(npcs))
        loadings = pca.fit_transform(tp_ref_vals)
        X = pca.components_.reshape(int(npcs), -1)
        varexpl = pca.explained_variance_ratio_
    else:
        loadings = []
        X = []
        varexpl = []
    
    return offset, pd.DataFrame(X, columns=tp_ref.columns, index=['PC_'+str(i) for i in range(int(npcs))]), pd.DataFrame(loadings, index=tp_ref.index, columns=['PC_'+str(i) for i in range(int(npcs))]), varexpl

def run_all_pcas(df_ref, npcs=1, timepoint_values=None, classes=["mesodermal", "neural"]):
    print(npcs)
    offset_dict = {}
    X_dict = {}
    Xvals_dict = {}
    loadings_dict = {}
    varexpl_dict = {}
    if type(npcs) is dict:
        max_npcs = int(max(npcs.values()))
    for cl in classes:
        if type(npcs) is dict:
            n = int(npcs[cl])
        else:
            n = int(npcs)
        offset_dict[cl], X_dict[cl], loadings_dict[cl], varexpl_dict[cl] = run_pca(df_ref, npcs=n, class_values=[cl], timepoint_values=timepoint_values) 
        Xvals_dict[cl] = X_dict[cl].values
        if type(npcs) is dict:
            Xvals_dict[cl] = np.vstack([Xvals_dict[cl], np.zeros([max_npcs - n, Xvals_dict[cl].shape[1]])])

    offset = np.stack(list(offset_dict.values()), axis=-1)
    X = np.stack(list(Xvals_dict.values()), axis=-1)
    if X.size == 0:
        X = None
    
    return offset, X, offset_dict, X_dict, loadings_dict, varexpl_dict

def run_model(toc,
              scSigs,
              offsets,
              #init_log_exposure=-10, # in code say set to -10000, not using this
              init_timings=0,
              learn_rate=0.01,
              poll_interval=100,
              max_it= 1e7,
              log_likelihood_tolerance=1e-6,
              sparsity_tolerance=1e-4
             ) :
    
    ####################
    # Define the model #
    ####################
    p,n = toc.shape
    if scSigs is not None:
        pcs, _, s = scSigs.shape
    else:
        _, s = offsets.shape
    #########
    # Inputs
    if scSigs is not None:
        S = tf.constant(scSigs.astype('float32'),name='fixed_signals')
    b = tf.constant(offsets.astype('float32'),name='fixed_offsets')
    k = tf.constant(toc.astype('float32'),name='obs_counts')
    ############
    # Exposures
    #These are the things we actually train
    #Initialise to -10000, which converts to 0.  This value will not be changed when fitting null
    zinit = np.random.uniform(0,1, size=[s,n])
    zinit = np.log((zinit/np.sum(zinit,axis=0))*np.sum(toc,axis=0))
    z = tf.Variable(zinit.astype('float32'),name='exposures')
    if scSigs is not None:
        x = tf.Variable(tf.zeros([pcs,n,s])+init_timings,name='timings')
    #Positive exposures
    E = tf.exp(z,name='positive_exposures')
    ###################
    # Predicted Counts
    if scSigs is not None:
        y = tf.matmul((tf.exp(tf.reduce_sum(x*S, axis=0)+b)),E,name='pred_counts')
    else:
        y = tf.matmul((tf.exp(b)),E,name='pred_counts')
    ##########################
    # Poisson log-likelihood
    #Variable part of Poisson log-likelihood
    Dij = k*tf.math.log(y)- y #obs*log(pred) - pred, neg of NLL because this is just LL
    #Constant part
    D0 = tf.math.lgamma(k+1)
    #Add gene weights and sum to get negative log-likelihood
    LL = tf.transpose(a=D0-Dij)
    tLL = tf.reduce_sum(input_tensor=LL,name='NLL')
    #############################
    # Define objective function
    #Count of non-zeroish co-efficients
    #cCnt = tf.reduce_sum(input_tensor=tf.sigmoid(zz))
    cCnt = tf.reduce_sum(input_tensor=tf.sigmoid(z))
    # The final penalised, adjust NLL
    O = tLL
    ############
    # Optimiser
    opt_E = tf.compat.v1.train.AdamOptimizer(learn_rate)
    
    ############
    # Main fit #
    ############
    print('''
    ######################
    # Fitting main Model #
    ######################
    ''')
    #Define optimisers
    if scSigs is not None:
        toLearn = [z,x]
    else:
        toLearn = [z]
    learners = opt_E.minimize(O,var_list=toLearn,name='learn_exposures')
    #Initialise
    sess = tf.compat.v1.Session()
    init = tf.compat.v1.global_variables_initializer()
    sess.run(init)
    #Record initial values
    last = sess.run(O)
    lastNonZero = sess.run(cCnt)
    #Record the movements
    nll = np.zeros(int(ceil(max_it/poll_interval)))
    nsigs = np.zeros(int(ceil(max_it/poll_interval)))
    i=0
    while True:
      #Take the exposure step
      sess.run([learners])
      #Every now and then, record our progress and check if we've converged
      if i%poll_interval == 0:
        ii = i//poll_interval
        #Record object function and number of non-zero exposures
        nll[ii] = sess.run(O)
        nsigs[ii] = sess.run(cCnt)
        #Record how much we've changed since we last checked
        diff = (last-nll[ii])
        last = nll[ii]
        diffCnts = (lastNonZero - nsigs[ii])
        lastNonZero = nsigs[ii]
        #Calculate per-element summaries
        sigsPerSample = lastNonZero/n
        llPerEntry = nll[ii]/n/p
        #And the average intercept
        #avgInt = np.mean(np.exp(sess.run(int0)))
        #The average coverage relative to the observed
        avgCov = np.mean(sess.run(y).sum(axis=0)/sess.run(k).sum(axis=0))
        print("[%s] step %d, training O=%g, cnt=%g,dO=%g,dCnt=%g,nSigsAvg=%g/%d,avgNLL=%g,avgCov=%g" %(datetime.datetime.now(),i,last,lastNonZero,diff,diffCnts,sigsPerSample,s+1,llPerEntry,avgCov))
        #Test if we should terminate
        if np.isnan(last):
            break
        if diff<0 and diff/last < log_likelihood_tolerance  and diffCnts/lastNonZero < sparsity_tolerance:
            break
      i = i+1
      if i>max_it:
        break
    
    E_out = sess.run(E)
    if scSigs is not None:
        x_out = sess.run(x)
    else:
        x_out = None
    y_pred = sess.run(y)
    
    return E_out, x_out, y_pred, nll[nll > 0]