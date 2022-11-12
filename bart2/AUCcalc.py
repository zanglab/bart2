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
import math
# import multiprocessing # multiprocessing on dealing with TF datasets

def get_tf_file_data(tf_json):
    starttime = time.time()
    with open(tf_json, 'r') as fm:
        tf_file_map = json.load(fm)
    endtime = time.time()
    sys.stdout.write("Loading TR file mapping list: {} seconds \n".format(endtime-starttime))
    return tf_file_map
    
# 84M memory taken
def get_matrix_data(overlap_json):
    starttime = time.time()
    with open(overlap_json, 'r') as fm:
        matrix_data = json.load(fm)
    endtime = time.time()
    sys.stdout.write("Loading TR matrix file: {} seconds \n ".format(endtime-starttime))
    return matrix_data

def get_binding_count(binding_count_json):
    starttime = time.time()
    with open(binding_count_json, 'r') as fm:
        bc_data = json.load(fm)
    endtime = time.time()
    sys.stdout.write("Loading TR matrix file: {} seconds \n ".format(endtime-starttime))
    return bc_data

#find estimation bins in the (left,right) interval
def find_est_bin(left, right, bin_size):
    if right//bin_size == left//bin_size:
        bins_insert = [(left,right)]
    else:
        #number of est boundries
        nbds=(right+1)//bin_size-(left-1)//bin_size
        #list of boundries (left-closed, right-open)
        bds=[left]
        for j in range(0,nbds):
            bds.append(((left-1)//bin_size+1+j)*bin_size)
        bds.append(right+1)
        #remove duplicated boundries if there are (could happen at left or right end)
        bds=sorted(list(set(bds)))
        #return bins to insert
        bins_insert=[]
        for j in range(0,len(bds)-1):
            bins_insert.append((bds[j],bds[j+1]-1))
    return(bins_insert)

def cal_bin_auc(bin_list, baseA, matrix_data, positions):
    start = bin_list[0]
    end = bin_list[1]
    y_addition = {}
    baseB = {}
    bin_auc = {}
    #initiate dicts
    for i in range(1, len(baseA)+1):
        y_addition[i] = 0
        bin_auc[i] = 0
        baseB[i] = 0
    #caculate addition on y-axis in the given bin for each TF
    for i in range(start,end+1):
        tf_occur = matrix_data[str(positions[i])].strip().split()
        for j in {int(t) for t in tf_occur }:
            y_addition[j] += 1
    #caculate baseB and AUC
    for i in range(1, len(baseA)+1):
        baseB[i] = baseA[i] + y_addition[i]
        bin_auc[i] = (baseA[i]+baseB[i])*(end-start+1-y_addition[i])/2
    return((bin_auc,baseB))


def cal_auc_for_all_tfs(args, positions, active_pos_count, tied_list, matrix_data, bc_data, tf_file_len):
    starttime0 = time.time()
    if args.subcommand_name == 'region':
        print('region mode - active positions:'+str(active_pos_count))
        active_positions = positions[0:active_pos_count]
        
        print("Determining bins for estimation...")
        starttime = time.time()
        #0-initialize bin intervals for estimation
        est_bin_size=100
        
        #merge est bins and tied bins
        #1) remove all tied bins fully contained by any est bin
        tied_list_f1=[]
        for tb in tied_list:
            if tb[0]//est_bin_size != tb[1]//est_bin_size:
                tied_list_f1.append(tb)

        #2) go through tied bins and estimate bins to merge them together
        all_bins = []

        if len(tied_list_f1) == 0:
            all_bins = find_est_bin(0, len(positions)-1, est_bin_size)
        else:
            for i in range(0,(len(tied_list_f1)-1)):
                if tied_list_f1[i][1]>=tied_list_f1[i+1][0]:
                    sys.exit("interval division error")
                elif tied_list_f1[i][1]+1==tied_list_f1[i+1][0]:
                    #neighouring tied bins connect
                    all_bins.append(tied_list_f1[i])
                else:
                    bins_insert = find_est_bin(tied_list_f1[i][1]+1, tied_list_f1[i+1][0]-1, est_bin_size)
                    all_bins.append(tied_list_f1[i])
                    all_bins = all_bins + bins_insert
            all_bins.append(tied_list_f1[-1])
            if tied_list_f1[0][0]!=0:
                head_bins = find_est_bin(0, (tied_list_f1[0][0]-1), est_bin_size)
                all_bins = head_bins + all_bins
            if tied_list_f1[-1][1]!=len(positions)-1:
                tail_bins = find_est_bin(tied_list_f1[-1][1], len(positions)-1, est_bin_size)
                all_bins = all_bins + tail_bins

        #check bin division
        for i in range(0,(len(all_bins)-1)):
            if all_bins[i][1]+1 != all_bins[i+1][0]:
                sys.exit("interval division error")
        if all_bins[0][0] != 0 or all_bins[-1][1] != len(positions)-1:
            sys.exit("interval division error")

        endtime = time.time()
        sys.stdout.write("{} seconds \n".format(endtime-starttime))       
        print("Divided all UDHS into "+str(len(all_bins))+" bins")

        #area calculation
        #calculate area units for each bin
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

        #multiply area units with unit size
        tf_auc = {}
        for i in range(1, tf_file_len+1):
            tf_auc[i] = tf_auc_unit[i]/(baseA[i]*(len(positions)-baseA[i]))

        endtime0 = time.time()
        sys.stdout.write("cal_auc_for_all_tfs takes {} seconds \n".format(endtime0-starttime0))

        return(tf_auc)


        #step1: calculate the area units of active part
        starttime = time.time()
        print("region mode: step1...")
        tf_auc_ac = {}
        height = {}#final heights are for step2
        tf_set = set(range(1, tf_file_len+1))
        #initiation
        for i in range(1, tf_file_len+1):
            height[i] = 0
            tf_auc_ac[i] = 0
        #calculate
        loop_count=0
        for i in active_positions:
            tf_occur = set(matrix_data[str(i)].strip().split())
            tf_occur = {int(t) for t in tf_occur }
            for j in tf_set-tf_occur:
                tf_auc_ac[j]+=height[j]
            for j in tf_occur:
                height[j]+=1
            loop_count+=1
            if loop_count%10000==0:
                print(loop_count)
        endtime = time.time()
        sys.stdout.write("{} seconds \n".format(endtime-starttime))

        #step2: calculate the area units of incative part (trapezoid)
        #base1:height of step 1 base2:binding count height: 2.7M-binding count-x axis moves in step in
        print("region mode: step2...")
        tf_auc_inac = {}
        for i in range(1, tf_file_len+1):
            tf_auc_inac[i]=(height[i]+bc_data[str(i)]) * (len(positions)-bc_data[str(i)]-(active_pos_count-height[i])) / 2
        endtime = time.time()

        #step3: add two parts and divided by total count of units
        print("region mode: step3...")
        tf_auc = {}
        for i in range(1, tf_file_len+1):
            tf_auc[i] = (tf_auc_ac[i] + tf_auc_inac[i]) / (bc_data[str(i)] * (len(positions)-bc_data[str(i)]))
        endtime = time.time()

        endtime0 = time.time()
        sys.stdout.write("cal_auc_for_all_tfs takes {} seconds \n".format(endtime0-starttime0))
        
    else:
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

def cal_auc(args, positions, active_pos, tied_list):
    tf_json = args.tffile
    overlap_json = args.tfoverlap
    normfile = args.normfile
    binding_count_json = args.bc
 
    # in case position index is not 'int' type
    positions = [int(i) for i in positions]
    if len(positions) == 0:
        sys.stderr.write('Input file might not with right format!\n')
        sys.exit(1)
    sys.stdout.write("Calculating ROC-AUC values for all transcriptional regulators:\n\n")

    tf_dict = get_tf_file_data(tf_json)
    overlap_dict = get_matrix_data(overlap_json)
    bc_dict = get_binding_count(binding_count_json)
    tf_auc = cal_auc_for_all_tfs(args, positions, active_pos, tied_list, overlap_dict, bc_dict, len(tf_dict))

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
