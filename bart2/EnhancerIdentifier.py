import os,sys,argparse
import tables
import time,re
import numpy as np

from bart2 import RPRegress

# read in H3K27ac sample IDs
def read_sample_list(fname):
    fp = open(fname)
    sl,scores = [],[]
    for line in fp:
        line  = line.strip().split('\t')
        #print(line)
        if re.match('[0-9]',line[0]):
            sl += [line[0]]
            scores += [float(line[1])]
    fp.close()
    return sl,scores

def median_norm(X):
    col_median = np.median(X,0)
    for i in range(X.shape[1]):
        X[:,i] = X[:,i] - col_median[i]
    return X

def main(args):
    '''generate cis-regulatory profile from selected H3K27ac samples

    Import arguments:
    argparser arguments

    Output file:
    Enhancer profile based on UDHS or UAC.
    '''
    samplefile = args.samplefile
    output_name = args.name
    H3K27ac_hdf5 = args.k27achdf5

    sample_names,sample_weights = read_sample_list(samplefile)
    sample_names = [i for i in sample_names if re.match('[0-9]',i)] # incase 'start/end/chr' in sample names
    sample_weights = np.array(sample_weights)
    
    DHS_sample_names = [ elem+'_Strength' for elem in sample_names ]
    print("Selected H3K27ac samples...")
    print(DHS_sample_names)
    sys.stdout.flush()

    # read data from RPKM H3K27ac hdf5 file
    udhs_h5file = tables.open_file(H3K27ac_hdf5, driver="H5FD_CORE")    
    chrom = RPRegress.read_hdf5(udhs_h5file, "chrom")
    start = RPRegress.read_hdf5(udhs_h5file, "start")
    end = RPRegress.read_hdf5(udhs_h5file, "end")
    ID = RPRegress.read_hdf5(udhs_h5file, "ID")
    DHS = None
    for dhs_samplename in DHS_sample_names:
        dhs = RPRegress.read_hdf5(udhs_h5file, dhs_samplename)
        dhs = np.array(dhs)
        dhs = dhs.transpose()
        if DHS is None:
            DHS = dhs
        else:
            DHS =  np.vstack((DHS,dhs))
    udhs_h5file.close()
    DHS = DHS.transpose()

    DHS = RPRegress.sqrttransform(DHS)
    T = np.dot(DHS,sample_weights) # coef * H3K27ac RPKM signals

    out_res = []
    for i,elem in enumerate(T):
        out_res.append((chrom[i].decode('utf-8'),start[i],end[i],str(i+1),elem))
    sorted_out_res = sorted(out_res, key=lambda x: float(x[4]),reverse=True)
    
    fpo = open(output_name+'_enhancer_prediction_lasso.txt','w')
    print('%s\t%s\t%s\t%s\t%s' % ("chromosom",'start','end','UDHSID',"Score"), file=fpo)
    for line in sorted_out_res:
        print('%s\t%d\t%d\t%s\t%3.2f' % (line[0],line[1],line[2],line[3],line[4]), file=fpo)
    fpo.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Predict enhancer elements from gene list and H3K27ac in union DNase-seq peaks.""")
    parser.add_argument( '-s', dest='samplefile', type=str, required=True, help='File that lists informative H3K27ac samples. One sample ID per line.' )
    parser.add_argument( '--k27ac', dest='k27achdf5', type=str, required=True, help='Path for the hdf5 format H3K27ac reads counts in UDHS regions.' )
    parser.add_argument( '-n', dest='name', type=str, required=True, help='Name of study, for output file naming.' )
    # parser.add_argument( '-t', dest='tffile', default=None, required=False, help='Indicators of TF binding at each UDHS sites 0 or 1. For performance evaluation.' )
    
    args = parser.parse_args()
    main(args)
