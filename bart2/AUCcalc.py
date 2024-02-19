# Time-stamp: <2017-08-10>
'''Module for calculating ROC-AUC values for all TF datasets

Copyright (c) 2017, 2018 Chongzhi Zang, Zhenjia Wang <zhenjia@virginia.edu>

This code is free software; you can redistribute it and/or modify it 
under the terms of the BSD License.

@status: release candidate
@version: $Id$
@author: Chongzhi Zang, Zhenjia Wang
@contact: zhenjia@virginia.edu

'''

import os
import sys
import time
import json
import math
from bart2 import WriteFile
# import multiprocessing # multiprocessing on dealing with TF datasets


def get_tf_file_data(tf_json):
    starttime = time.time()
    with open(tf_json, 'r') as fm:
        tf_file_map = json.load(fm)
    endtime = time.time()
    sys.stdout.write(
        "Loading TR file mapping list: {} seconds \n".format(endtime-starttime))
    return tf_file_map

# 84M memory taken


def get_matrix_data(overlap_json):
    starttime = time.time()
    with open(overlap_json, 'r') as fm:
        matrix_data = json.load(fm)
    endtime = time.time()
    sys.stdout.write(
        "Loading TR matrix file: {} seconds \n ".format(endtime-starttime))
    return matrix_data

# find estimation bins within the (left,right) interval (closed interval)


def find_est_bin(left, right, bin_size):
    if right//bin_size == left//bin_size:
        bins_insert = [(left, right)]
    else:
        # number of est boundries
        nbds = (right+1)//bin_size-(left-1)//bin_size
        # list of boundries (left-closed, right-open)
        bds = [left]
        for j in range(0, nbds):
            bds.append(((left-1)//bin_size+1+j)*bin_size)
        bds.append(right+1)
        # remove duplicated boundries if there are (could happen at left or right end)
        bds = sorted(list(set(bds)))
        # return bins to insert
        bins_insert = []
        for j in range(0, len(bds)-1):
            bins_insert.append((bds[j], bds[j+1]-1))
    return (bins_insert)


def cal_bin_auc(bin_list, baseA, matrix_data, positions):
    start = bin_list[0]
    end = bin_list[1]
    y_addition = {}
    baseB = {}
    bin_auc = {}
    # initiate dicts
    for i in range(1, len(baseA)+1):
        y_addition[i] = 0
        bin_auc[i] = 0
        baseB[i] = 0
    # caculate addition on y-axis in the given bin for each TF
    for i in range(start, end+1):
        tf_occur = matrix_data[str(positions[i])].strip().split()
        for j in {int(t) for t in tf_occur}:
            y_addition[j] += 1
    # caculate baseB and AUC
    for i in range(1, len(baseA)+1):
        baseB[i] = baseA[i] + y_addition[i]
        bin_auc[i] = (baseA[i]+baseB[i])*(end-start+1-y_addition[i])/2
    return ((bin_auc, baseB))


def cal_auc_for_all_tfs(args, positions, tied_list, matrix_data, tf_file_len):
    starttime0 = time.time()

    print("Determining bins for estimation...")
    starttime = time.time()
    # bin size for estimation
    if args.binsize is None:
        if args.subcommand_name == 'geneset':
            est_bin_size = 1000
        elif args.subcommand_name == 'profile':
            est_bin_size = 1000
        elif args.subcommand_name == 'region':
            est_bin_size = 50
    else:
        est_bin_size = args.binsize

    # merge est bins and tied bins
    # 1) remove all tied bins fully contained by any est bin
    tied_list_f1 = []
    for tb in tied_list:
        if tb[0]//est_bin_size != tb[1]//est_bin_size:
            tied_list_f1.append(tb)

    # 2) go through tied bins and estimate bins to merge them together
    all_bins = []

    if len(tied_list_f1) == 0:
        all_bins = find_est_bin(0, len(positions)-1, est_bin_size)
    else:
        for i in range(0, (len(tied_list_f1)-1)):
            if tied_list_f1[i][1] >= tied_list_f1[i+1][0]:
                sys.exit("interval division error")
            elif tied_list_f1[i][1]+1 == tied_list_f1[i+1][0]:
                # neighouring tied bins connect
                all_bins.append(tied_list_f1[i])
            else:
                bins_insert = find_est_bin(
                    tied_list_f1[i][1]+1, tied_list_f1[i+1][0]-1, est_bin_size)
                all_bins.append(tied_list_f1[i])
                all_bins = all_bins + bins_insert
        all_bins.append(tied_list_f1[-1])
        if tied_list_f1[0][0] != 0:
            head_bins = find_est_bin(0, (tied_list_f1[0][0]-1), est_bin_size)
            all_bins = head_bins + all_bins
        if tied_list_f1[-1][1] != len(positions)-1:
            tail_bins = find_est_bin(
                tied_list_f1[-1][1]+1, len(positions)-1, est_bin_size)
            all_bins = all_bins + tail_bins

    # check bin division
    for i in range(0, (len(all_bins)-1)):
        if all_bins[i][1]+1 != all_bins[i+1][0]:
            print(all_bins[i])
            print(all_bins[i+1])
            sys.exit("interval division error")
    if all_bins[0][0] != 0 or all_bins[-1][1] != len(positions)-1:
        sys.exit("interval division error")

    endtime = time.time()
    sys.stdout.write("{} seconds \n".format(endtime-starttime))
    print("Divided all UDHS into "+str(len(all_bins))+" bins")

    # area calculation
    # calculate area units for each bin
    baseA = {}
    tf_auc_unit = {}
    for i in range(1, tf_file_len+1):
        tf_auc_unit[i] = 0
        baseA[i] = 0

    for i in range(0, len(all_bins)):
        bin_area_res = cal_bin_auc(all_bins[i], baseA, matrix_data, positions)
        for j in range(1, tf_file_len+1):
            tf_auc_unit[j] += bin_area_res[0][j]
            baseA[j] = bin_area_res[1][j]

    # multiply area units with unit size
    tf_auc = {}
    for i in range(1, tf_file_len+1):
        tf_auc[i] = tf_auc_unit[i]/(baseA[i]*(len(positions)-baseA[i]))

    endtime0 = time.time()
    sys.stdout.write(
        "cal_auc_for_all_tfs takes {} seconds \n".format(endtime0-starttime0))

    return (tf_auc)


def get_position_list(enhancerfile):
    '''
    Get the ID list of DHS, according to the decreasingly sorted scores in MARGE enhancer profile
    '''
    fin = open(enhancerfile, 'rb')
    line = fin.readline()
    score = {}
    while line:
        line = line.strip().split()
        try:
            score[line[-2]] = float(line[-1])
        except:
            pass
        line = fin.readline()
    fin.close()
    return sorted(score.keys(), key=score.get, reverse=True)


def cal_auc(args, positions, tied_list):
    tf_json = args.tffile
    overlap_json = args.tfoverlap
    normfile = args.normfile

    # in case position index is not 'int' type
    positions = [int(i) for i in positions]
    if len(positions) == 0:
        sys.stderr.write('Input file might not with right format!\n')
        sys.exit(1)
    sys.stdout.write(
        "Calculating ROC-AUC values for all transcriptional regulators:\n\n")

    tf_dict = get_tf_file_data(tf_json)
    overlap_dict = get_matrix_data(overlap_json)
    tf_auc = cal_auc_for_all_tfs(
        args, positions, tied_list, overlap_dict, len(tf_dict))

    return tf_auc, tf_dict


if __name__ == '__main__':
    margefile = '/nv/vol190/zanglab/wm9tr/software/BART-v1.0.1-py3-full/BART/hg38_library/hg38_test_data/hg38_AR_up_genes_enhancer_prediction.txt'
    tf_json = '/nv/vol190/zanglab/wm9tr/test/test_bart_storage/human_file_str_format.json'
    overlap_json = '/nv/vol190/zanglab/wm9tr/test/test_bart_storage/human_matrix_str_format.json'

    # tf_json = '/nv/vol190/zanglab/wm9tr/test/test_bart_storage/human_UAC_file_str_format.json'
    # overlap_json = '/nv/vol190/zanglab/wm9tr/test/test_bart_storage/human_UAC_matrix_str_format.json'

    output_name = 'test_bart/AR_aupr'
    norm_file = '/nv/vol190/zanglab/wm9tr/software/BART-v1.0.1-py3-full/BART/hg38_library/hg38_MSigDB.dat'

    cal_auc(args, positions)
