# Time-stamp: <2017-08-10>
'''Module for calculating the Wilcoxon-score and p value for each unique TF 

Copyright (c) 2017, 2018 Chongzhi Zang, Zhenjia Wang <zhenjia@virginia.edu>

This code is free software; you can redistribute it and/or modify it 
under the terms of the BSD License.

@status: release candidate
@version: $Id$
@author: Chongzhi Zang, Zhenjia Wang
@contact: zhenjia@virginia.edu

'''
import os,sys,re
import pandas as pd
import numpy as np
import scipy
from scipy import stats

def factorial(n):
    value = 1.0
    while n>1:
        value*=n
        n-=1
    return value

def logfac(n):
    if n<20:
        return np.log(factorial(n))
    else:
        return n*np.log(n)-n+(np.log(n*(1+4*n*(1+2*n)))/6.0)+(np.log(np.pi))/2.0

def irwin_hall_cdf(x,n):
    # pval = returned_value for down regulated
    # pval = 1 - returned_value for up regulated
    value,k = 0,0
    while k<=np.floor(x):
        value +=(-1)**k*(scipy.special.binom(n,k))*(x-k)**n
        k+=1
    return value/(np.exp(logfac(n)))


# def stat_test(AUCs,args): 
def stat_test(AUCs, tf_dict, statfile, normfile): 
    # read AUCs according to TF type
    sys.stdout.write('Statistical tests start.\n')
    tfs = {}
    sam1 = []
    for tf_key in AUCs.keys():
        tf = tf_dict[str(tf_key)].split('_')[0]
        auc = AUCs[tf_key]
        sam1.append(auc)
        if tf not in tfs:
            tfs[tf] = [auc]
        else:
            tfs[tf].append(auc)
    
    cols = ['score','pvalue','max_auc','zscore','rank_score','rank_zscore','rank_pvalue','rank_auc','rank_avg_z_p','rank_avg_z_p_a','rank_avg_z_p_a_irwinhall_pvalue']
    stat = pd.DataFrame(index = [tf for tf in tfs],columns = cols)
    for tf in tfs.keys():
        if len(tfs[tf])>0: # filter the tf with few samples
            stat_test = stats.ranksums(tfs[tf],sam1)
            stat.loc[tf]['score'] = stat_test[0]
            # one-sided test
            stat.loc[tf]['pvalue'] = stat_test[1]*0.5 if stat_test[0]>0 else 1-stat_test[1]*0.5

    tf_stats = pd.read_csv(normfile, sep='\t', index_col=0)
    # cal the normalized stat-score 
    sys.stdout.write('Do standardization...')
    for i in stat.index:
        #stat[i].append((stat[i][0]-tf_stats.loc[i,'mean'])/tf_stats.loc[i,'std']) #[2] for Z-Score
        try:
            stat.loc[i]['zscore'] = (stat.loc[i]['score']-tf_stats.loc[i,'mean'])/tf_stats.loc[i,'std']
            stat.loc[i]['max_auc'] = max(tfs[i])
        except KeyError:
            stat.loc[i]['zscore'] = 0.0
            stat.loc[i]['max_auc'] = 0.0

    # rank the list by the average rank of stat-score and z-score
    # rank of Wilcoxon Socre
    rs = 1
    for i in sorted(stat.index,key = lambda x: stat.loc[x]['score'],reverse=True): 
        stat.loc[i]['rank_score'] = rs #  rank of stat_score
        rs +=1

    # rank of Z-Score
    rz = 1
    for i in sorted(stat.index,key = lambda x: stat.loc[x]['zscore'],reverse=True):        
        stat.loc[i]['rank_zscore'] = rz # rank of z-score
        rz +=1 

    # rank of pvalue
    rp = 1
    for i in sorted(stat.index,key = lambda x: stat.loc[x]['pvalue'],reverse=False):        
        stat.loc[i]['rank_pvalue'] = rp #  rank of pvalue
        rp +=1

    ra = 1
    for i in sorted(stat.index,key = lambda x: stat.loc[x]['max_auc'],reverse=True):        
        stat.loc[i]['rank_auc'] = ra #  rank of pvalue
        ra +=1

    # rank of average
    for i in stat.index:
        stat.loc[i]['rank_avg_z_p'] = (stat.loc[i]['rank_zscore']+stat.loc[i]['rank_pvalue'])*0.5   # [6] for average of stat-score and z-score
        stat.loc[i]['rank_avg_z_p_a'] = (stat.loc[i]['rank_zscore']+stat.loc[i]['rank_pvalue']+stat.loc[i]['rank_auc'])*0.33/len(tfs.keys())   # [7] for average of three
        stat.loc[i]['rank_avg_z_p_a_irwinhall_pvalue'] = irwin_hall_cdf(3*stat.loc[i]['rank_avg_z_p_a'],3)

    with open(statfile,'w') as statout:
        statout.write('TF\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('statistic','pvalue','zscore','max_auc','re_rank','irwin_hall_pvalue'))
        for i in sorted(stat.index,key=lambda x: stat.loc[x]['rank_avg_z_p_a'],reverse=False):
            statout.write('{}\t{:.3f}\t{:.3e}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3e}\n'.format(i,stat.loc[i]['score'],stat.loc[i]['pvalue'],stat.loc[i]['zscore'],stat.loc[i]['max_auc'],stat.loc[i]['rank_avg_z_p_a'],stat.loc[i]['rank_avg_z_p_a_irwinhall_pvalue']))
    sys.stdout.write('--Standardization finished!\n--Ranked TFs saved in file: {}\n'.format(statfile))

    # plot figures of user defined TFs
    # if args.target:
    #     with open(args.target) as target_file:
    #         IDs = [re.split('[^a-zA-Z0-9]+',line)[0] for line in target_file.readlines()]
    #         for ID in IDs:
    #             stat_plot(stat,tfs,ID,args,'rank_avg_z_p_a')

    sys.stdout.write('Prediction done!\n')
