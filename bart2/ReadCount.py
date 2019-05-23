# Time-stamp: <2017-08-10>
'''
Copyright (c) 2017, 2018 Chongzhi Zang, Zhenjia Wang <zhenjia@virginia.edu>

This code is free software; you can redistribute it and/or modify it 
under the terms of the BSD License.

@status: release candidate
@version: $Id$
@author: Chongzhi Zang, Zhenjia Wang
@contact: zhenjia@virginia.edu

'''

import os,re
import bisect
from bart2 import IOparser

plus = re.compile('\+')
minus = re.compile('\-')
      

def is_list_sorted(mylist):
    '''
    Check if list is sorted
    '''        
    for i in range(len(mylist)-1):
        if mylist[i] > mylist[i+1]:
            return 0
    return 1   
        	

def get_read_positions(positions,regions,val,fragment_size):
    '''
    Return all the shifted positions of reads 
    '''
    # default fragment_size = 150
    if fragment_size >= 0:
        shift = int(round(fragment_size/2))	
        for chrom in regions.keys():
            if chrom not in positions:
                positions[chrom]=[]
            for inner in regions[chrom].keys():
                positions[chrom].extend([outer+shift*val for outer in regions[chrom][inner]])
    else:
        for chrom in regions.keys():
            if chrom not in positions:
                positions[chrom]=[]
            for inner in regions[chrom].keys():
                for outer in regions[chrom][inner]:
                    shift = int(abs(outer-inner)/2)	
                    positions[chrom].append(outer+shift*val)
    return positions   
	        	    
		
def get_count_on_DHS(start,end,positions):
    '''
    Count the tags/positions on DHS
    '''
    #if is_list_sorted(positions)==0:
    #    positions.sort()
    if start <= end:
        middle = int((start+end)*0.5)
        s = bisect.bisect_left(positions,middle-500)
        e = bisect.bisect_right(positions,middle+500)
        return e-s
    else:
        return 0


def read_count_on_DHS(args):
    '''
    Count the num of (unique) reads in user-input bed/bam file
    on each of the UDHS -- bart profile
    '''
    #specify the species
    # get the start-end regions of each read in each chrom (as key) and separate by strand

    regions1, regions2 = IOparser.get_tag_regions(args.species,args.format,args.infile)

    # get the counting position of each read(tag), chrom as key
    positions = get_read_positions({},regions1,1,args.fragmentsize)
    positions = get_read_positions(positions,regions2,-1,args.fragmentsize)
    total=0
    for chrom in positions:
        positions[chrom].sort() 
        total+=len(positions[chrom])
    # check if positions is NULL
    #if total ==0:
        #sys.stderr.write('Can not read the input bed/bam file!\n')
    # count reads on each DHS

    DHSfile = args.dhsfile#os.path.dirname(__file__)+os.sep+'data'+os.sep+'{}_unionDHS.bed'.format(args.species)#
    Dfile = open(DHSfile,'r')
    line = Dfile.readline()
    counting = {}
    while line:
        line = line.strip().split()
        #dhs = BED(line[0],line[1],line[2],line[3],line[4],line[5])
        chrom = line[0]
        start = int(line[1])
        end = int(line[2])
        DHS_id = line[3]
        if chrom in positions:
            #assert DHS_id not in counting
            nums = get_count_on_DHS(start,end,positions[chrom])
            counting[DHS_id]= round(nums*1000000000/(total*1000),3)
        else:
            counting[DHS_id] = 0.0
        line = Dfile.readline()
    Dfile.close()
    return counting
