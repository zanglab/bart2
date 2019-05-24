#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 21:11:46 2019

@author: yifan
"""

import os,re
from bart2 import IOparser


def get_score_on_DHS(UDHS_start,UDHS_end,region_score):
    for key_start in region_score:
        if key_start >= UDHS_start and key_start <= UDHS_end:
            for key_end in region_score[key_start]:
                return int(list(region_score[key_start][key_end])[0])
        else:
            for key_end in region_score[key_start]:
                if key_end >= UDHS_start and key_end <= UDHS_end:
                    return int(list(region_score[key_start][key_end])[0])
    return 0



def score_on_DHS(args):
    '''
    Assign the score on bed interval to the overlapping UDHS
    '''
    #specify the species
    # get the start-end and score of each bed interval, save in a dict of chr:{start:end:score}
    regions = IOparser.get_bed_score(args.infile,args.species)

    # sort each intervals by start position
    regions_sorted = {}
    for chrm in regions:
        if chrm not in regions_sorted:
            regions_sorted[chrm] = {}
        for start in sorted(regions[chrm]):
            regions_sorted[chrm][start] = regions[chrm][start]

    #read in the UDHS intervals
    DHSfile = args.dhsfile#os.path.dirname(__file__)+os.sep+'data'+os.sep+'{}_unionDHS.bed'.format(args.species)#
    Dfile = open(DHSfile,'r')
    line = Dfile.readline()
    #initiate the returen, a dict of UDHS_number:score_on_bed
    counting = {}
    while line:
        line = line.strip().split()
        #UDHS interval information
        chrom = line[0]
        start = int(line[1])
        end = int(line[2])
        DHS_id = line[3]
        if chrom in regions:
            #get the score on this DHS, 0 if not overlapping with inpit bed file
            nums = get_score_on_DHS(start,end,regions_sorted[chrom])
            #remove any entry in the sorted bed region that come before the current DHS (start position)
            bed_key = list(regions_sorted[chrom])
            for key_start in bed_key:
                if key_start < start:
                    del regions_sorted[chrom][key_start]
                else:
                    break
            #register the count on this DHS
            counting[DHS_id]= nums
        else:
            counting[DHS_id] = 0.0
        line = Dfile.readline()
    Dfile.close()
    return counting
