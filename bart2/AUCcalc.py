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

import os,sys,time
import json
# import multiprocessing # multiprocessing on dealing with TF datasets

def get_tf_file_data(tf_json):
    starttime = time.time()
    with open(tf_json, 'r') as fm:
        tf_file_map = json.load(fm)
    endtime = time.time()
    sys.stdout.write("Loading tf file mapping list: {} seconds \n".format(endtime-starttime))
    return tf_file_map
    
# 84M memory taken
def get_matrix_data(overlap_json):
    starttime = time.time()
    with open(overlap_json, 'r') as fm:
        matrix_data = json.load(fm)
    endtime = time.time()
    sys.stdout.write("Loading tf matrix file: {} seconds \n ".format(endtime-starttime))
    return matrix_data

def cal_auc_for_all_tfs(args, positions, matrix_data, tf_file_len):
    udhs_len = len(positions)
    groupsize = 10000   # 2.7M / 10000 = 272 steps
    groups = int(udhs_len/groupsize)

    sys.stdout.write('Total groups: {} \n'.format(groups))

    tf_auc = {} # each tf file auc
    tf_t = {} # each tf file, how many 1s
    for i in range(1, tf_file_len+1): # initiate the axis for each tf file
        tf_t[i] = []
        tf_auc[i] = 0.0
    
    #TODO: is it needed to multiprocess the data? tf_t for each group instead of a tf_t_accumulation?
    # if args.processes:
    #     sys.stdout.write('--Number of cores will be used: {}\n'.format(args.processes))
    #     pool = multiprocessing.Pool(processes=args.processes)

    sys.stdout.write("Parsing data by group size: {} \n".format(groupsize))
    sys.stdout.flush()
    for group_id in range(groups+1):
        tf_t_cnt = {}
        for i in range(1, tf_file_len+1):
            tf_t_cnt[i] = 0

        for i in range(group_id*groupsize, min((group_id+1)*groupsize, udhs_len)):
            for tf in matrix_data[str(positions[i])].strip().split(): # count how many 1s in this group
                tf_t_cnt[int(tf)] += 1

        for tf, t_cnt in tf_t_cnt.items():
            if tf_t[tf]:
                tf_t[tf].append(tf_t[tf][-1]+t_cnt)
            else:
                tf_t[tf].append(t_cnt)

    for key, tf_t_cnt in tf_t.items():
        for i in range(len(tf_t_cnt)):
            cur_x = (min(groupsize*(i+1), udhs_len)-tf_t_cnt[i])/(udhs_len-tf_t_cnt[-1])
            cur_y = tf_t_cnt[i]/tf_t_cnt[-1]
            if i == 0:
                width = cur_x
                height = cur_y/2
                tf_auc[key] += height*width
            else:
                width = cur_x - (min(groupsize*(i), udhs_len)-tf_t_cnt[i-1])/(udhs_len-tf_t_cnt[-1])
                height = (cur_y + tf_t_cnt[i-1]/tf_t_cnt[-1])/2
                tf_auc[key] += height*width

    # TODO: ========below is for AUPR =====
    # for key, tf_t_value in each_tf_t.items():
    #     for i in range(len(tf_t_value)):
    #         cur_x = tf_t_value[i]/tf_t_value[-1]
    #         cur_y = tf_t_value[i]/min(groupsize*(i+1), udhs_len)
    #         if i == 0:
    #             width = cur_x
    #             height = cur_y
    #             tf_auc[key] += height*width
    #         else:
    #             width = cur_x - tf_t_value[i-1]/tf_t_value[-1]
    #             height = cur_y
    #             tf_auc[key] += height*width
    # ==========  end modification =========

    return tf_auc

def get_position_list(enhancerfile):
    '''
    Get the ID list of DHS, according to the decreasingly sorted scores in MARGE enhancer profile
    ''' 
    fin = open(enhancerfile,'rb')
    line = fin.readline()  
    score = {}
    while line:
        line = line.strip().split()
        try:
            score[line[-2]]=float(line[-1])
        except:
            pass
        line = fin.readline()
    fin.close()
    return sorted(score.keys(),key=score.get,reverse=True)

def cal_auc(args, positions):
    tf_json = args.tffile
    overlap_json = args.tfoverlap
    normfile = args.normfile
 
    # in case position index is not 'int' type
    positions = [int(i) for i in positions]
    if len(positions) == 0:
        sys.stderr.write('Input file might not with right format!\n')
        sys.exit(1)
    sys.stdout.write("Calculating ROC-AUC values for all transcription factors:\n\n")

    tf_dict = get_tf_file_data(tf_json)
    overlap_dict = get_matrix_data(overlap_json)
    tf_auc = cal_auc_for_all_tfs(args, positions, overlap_dict, len(tf_dict))

    auc_file = args.ofilename + '_auc.txt'
    # output file of AUC-ROC values for all TFs
    with open(auc_file, 'w') as aucf:
        for key in sorted(tf_auc.keys(),key=tf_auc.get,reverse=True):
            tf_file_name = tf_dict[str(key)]
            aucf.write('{}\tAUC = {:.3f}\n'.format(tf_file_name, tf_auc[key]))
    sys.stdout.write('\n--ROC-AUC calculation finished!\n--Results saved in file: {}\n'.format(auc_file))
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
