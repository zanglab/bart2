import argparse,math,os,sys,tables
import numpy as np
from operator import itemgetter
from sklearn import linear_model, model_selection, metrics
import random

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


UNCHANGED,DOWN,UP,TARGET = 0,1,2,3
statdict = { UNCHANGED:'.', DOWN:'DOWN', UP:'UP', TARGET:'TARGET' }

class dataset(object):
    def __init__(self):
        self.index = None
        self.info  = None # include chrom and TSS and gene symbol
        self.x     = None
        self.y     = None
        self.bwdir = None 
        self.bwfiles = None
        self.rpfiles = None 

def gene_sym(symfile):
    """
    One representative of each gene symbol is chosen. 

    """
    # TODO: change the selection from the first to the longest gene body

    # return {"refseq":["symbol","chr"]}
    fp = open(symfile)
    symdict = {} # {"gene_symbol": ["refseq_ID", "chr"]}  the first refseq ID in the TSS file

    for line in fp:
        f = line.strip().split('\t')
        g = f[3]
        IDs = g.split(':')
        if IDs[1] not in symdict:
            symdict[IDs[1]] = [IDs[0], f[0]]
 
    rsymdict = {}
    for elem in symdict:
        rsymdict[symdict[elem][0]] = [elem,symdict[elem][1]]
    fp.close()

    return rsymdict

def logtransform(x):
    xt = np.array(x)
    pcount = 1
    xt += pcount 
    med = np.median( xt,0 )
    x = np.log2(xt) - np.log2(med)
    return x

def sqrttransform(x):
    xt = np.array(x)
    med = np.median( xt,0 )
    x = np.sqrt(xt) - np.sqrt(med)
    return x

def read_hdf5( h5file, sample_name):
    """ Apply motif stat function to all data in motif_file_name. 
    Data in numpy array val corresponds to idx entries. If idx if None all entries are used.""" 
    a = h5file.get_node("/", sample_name )
    m = a.read()
    return m

def getSampleNames_hdf5(h5file):
    samplenames = []
    for array in h5file.walk_nodes("/","Array"):
        if array.name not in samplenames:
            samplenames.append(array.name)
        else:
            continue
    #samplenames = samplenames[340:]
    return samplenames

def readregpotfiles(sym,genome,samplenames,h5file):

    # make list of regulatory potential file
    # read in regulatory potential from files in directory 

    index = None
    x = None
    nsamples = len(samplenames)

    refseqID = read_hdf5( h5file, 'refseqID' )
    print(refseqID)
    symID = read_hdf5(h5file, 'symbol')
    print(symID)
    chrom = read_hdf5(h5file, 'chr')
    start = read_hdf5(h5file, 'start')
    for k,name in enumerate(samplenames):
        if index == None:
            index = {}
            info = {}
            i = 0
            for j,geneid in enumerate(refseqID):
                geneid = geneid.decode("utf-8")  # gene symbol
                if geneid in sym:
                    symid = sym[geneid][0]
                    if symid not in index:
                        index[symid] = i
                        info[symid] = [chrom[j].decode("utf-8"),start[j]]
                        i += 1
            ngenes = len(index)
            x = np.zeros((ngenes,nsamples))
            print(np.shape(x))
            
        RP = read_hdf5( h5file, name )
        for i,geneid in enumerate(refseqID):
            geneid = geneid.decode("utf-8")
            if geneid in sym:
                symid = sym[geneid][0]
                rp = RP[i]
                try:
                    x[index[symid],k] = rp  ### float num was ignored here, e.g., 'chr', 'refseqID', 'start', 'symbol'
                except:
                    pass
    
    z         = dataset()  
    z.rpfiles = samplenames 
    z.x = x # x.shape = ngenes,nsamples  # x is RP not relative RP, change to relative RP
    print(np.median(x, axis=1))

    z.index = index # {symbol:'start position'}
    z.info  = info # {'symbol':[chr,start]}
    return z


def read_genelistOnly(sym, fname, index, exptype):
    
    status = np.zeros( len(index) )
    genenames = np.ndarray(shape=(len(index)),dtype=object)
    print(list(index.keys())[0:20])

    train_chroms = ['chr1','chr3','chr5','chr7','chr9','chr11','chr13','chr15','chr17','chr19','chr21']
    test_chroms = ['chr2','chr4','chr6','chr8','chr10','chr12','chr14','chr16','chr18','chr20','chr22']
    train_index = []
    test_index = []

    fp = open(fname).readlines()
    genes = [g.strip() for g in fp]

    allgenes = list(sym.keys())
    print(allgenes[0:20])

    for ag in allgenes:
        if exptype == 'Gene_Only':
            try:
                i = index[sym[ag][0]]
                if sym[ag][1] in train_chroms:
                    train_index.append(i)
                elif sym[ag][1] in test_chroms:
                    test_index.append(i)
                #print i
                if sym[ag][0] in genes:
                    #print sym[ag][0]
                    status[i] = TARGET
                else:
                    status[i] = UNCHANGED
                genenames[i] = ag
            except:
                continue
        else:
            try:
                i = index[sym[ag][0]]
                if sym[ag][1] in train_chroms:
                    train_index.append(i)
                elif sym[ag][1] in test_chroms:
                    test_index.append(i)
                else:
                    pass
                if ag in genes:
                    status[i] = TARGET
                else:
                    status[i] = UNCHANGED
                genenames[i] = ag
            except:
                continue

    print('file: %s\ttarget: %d\tunchanged: %d\n' % ( fname, sum( status == TARGET ), sum( status == UNCHANGED ) ))
    # print(genenames[0:20])
    return (genenames, status,train_index,test_index)  

def dataset_annotation(annotationf):
    # get the cell annotation for each datasetID
    inf = open(annotationf,'rU')
    ann = {}
    for line in inf:
        if line.startswith('datasetID'):
            pass
        else:
            line = line.strip().split('\t')
            ID = line[0] # dataset id -> GSM id
            info = [line[4],line[5],line[7]] # CellLineName, Tissue/Organ, DetailedTissue
            try:
                ann[ID] = info
            except:
                ann[ID] = 'NA'
    return ann

def lasso_test(x,y):
    # given x,y, return auc score
    LR_l1 = linear_model.LogisticRegression(penalty='l1', tol=0.01)
    LR_l1.fit(x,y)
    #np.mean( model_selection.cross_val_score( LR_l1, X_train, y_train, scoring='roc_auc', cv=5 ))
    yhat = LR_l1.predict_log_proba(x)
    fpr, tpr, thresholds = metrics.roc_curve(y, yhat[:,1], pos_label=1)
    auc = metrics.auc(fpr,tpr)
    selected_features = len([i for i in LR_l1.coef_[0] if i !=0])
    return auc,selected_features
    
def lasso_test_best_alpha(x,y,prename):
    # given x,y, return alpha used for adaptive lasso
    alphas = [i for i in np.logspace(-2,1,10)]
    alpha_cvs = []
    plt.figure()
    for alpha in alphas:
        LR_l1 = linear_model.LogisticRegression(penalty='l1', tol=0.01,fit_intercept=True,C=alpha);print(alpha)
        cvs_scores = model_selection.cross_val_score( LR_l1, x, y, scoring='roc_auc', cv=5 )
        alpha_cvs.append(cvs_scores)
        
        LR_l1.fit(x,y)
        yhat = LR_l1.predict_log_proba(x)
        fpr, tpr, thresholds = metrics.roc_curve(y, yhat[:,1], pos_label=1)
        auc = metrics.auc(fpr,tpr)
        selected_features = len([i for i in estimator.coef_[0] if i !=0])
        print(alpha,np.mean(cvs_scores),auc,selected_features)
        # plot the auc figs
        y_mean = np.mean(cvs_scores)
        y_err  = np.std(cvs_scores)
        plt.errorbar(alpha,y_mean,y_err,color='r',ecolor='grey',fmt='o',capsize=4)
    plt.ylim([0,1])
    plt.xscale('log')
    plt.savefig(prename+'_alpha_auc.png',bbox_inches='tight',pad_inches=0.1,transparent=True)
    plt.close()
    #alpha_cvs_mean = [i.mean() for i in alpha_cvs]
    #best_alpha = alphas[alpha_cvs_mean.index(max(alpha_cvs_mean))]
    return alphas,alpha_cvs


def best_alpha(x,y):
    # given x,y, return alpha used for adaptive lasso
    alphas = [i for i in np.logspace(-2,1,10)]
    #alphas = [0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2]
    alpha_cvs = []
    
    for alpha in alphas:
        LR_l1 = linear_model.LogisticRegression(penalty='l1', tol=0.01,fit_intercept=True,C=alpha)
        cvs_scores = model_selection.cross_val_score( LR_l1, x, y, scoring='roc_auc', cv=5 )
        alpha_cvs.append(cvs_scores) 
        print('  =best-alpha= ',alpha, '==mean-cvs==',np.mean(cvs_scores))
    alpha_cvs_mean = [i.mean() for i in alpha_cvs]
    best_alpha = alphas[alpha_cvs_mean.index(max(alpha_cvs_mean))]
   
    return best_alpha,max(alpha_cvs_mean)


def adaptive_lasso(x,y,samplefiles,name,maxsamples,ann,genenames):
    # test of adaptive lasso
    g = lambda w:np.sqrt(np.abs(w))
    gprime = lambda w: 1.0/(2.*np.sqrt(np.abs(w))+np.finfo(float).eps)
    n_samples,n_features = x.shape
    
    n_lasso_iterations = 10
    weights = np.ones(n_features)
    selected_features = n_features
    
    print('Run adaptive lasso for 10 rounds...')
    print('Round, Alpha, Features number ')
    sys.stdout.flush()
    for k in range(n_lasso_iterations):
        if selected_features >maxsamples: 
            alpha=0.02
        else:
            alpha=0.2
        X_w = x / weights
        #alpha,best_cvs = best_alpha(X_w,y) # TODO: if you need to select best alpha for each step later
        #alpha = 0.1

        # set fixed seed and default solver
        estimator = linear_model.LogisticRegression(penalty='l1', tol=0.01, fit_intercept=True, C=alpha, random_state=2019, solver="liblinear")
        estimator.fit(X_w,y)
        coef_ = estimator.coef_/weights
        weights = gprime(coef_)
        selected_features = len([i for i in coef_[0] if i !=0])
        print('{}, {}, {}'.format(k,alpha,selected_features))
        sys.stdout.flush()

    rand_idx = list(range(x.shape[0]))
    random.shuffle( rand_idx )
    # xt = np.multiply(x,coef_);print(xt.shape)
    xt,yt = X_w[rand_idx,:], y[rand_idx]
    cvs_scores = model_selection.cross_val_score(estimator ,xt,yt,  scoring='roc_auc', cv=5 )
    best_cvs = np.mean(cvs_scores)
    yhat = estimator.predict_log_proba(xt)
    fpr, tpr, thresholds = metrics.roc_curve(yt, yhat[:,1], pos_label=1)
    auc = metrics.auc(fpr,tpr)
    
#    print(k,'alpha',alpha)
#    print(k,'best_cvs',best_cvs)
#    print(k,'auc',auc)
#    print(k,'selected_features',selected_features)

    outf = open('{}_adaptive_lasso_Info.txt'.format(name),'w')
    for coef in sorted([ i for i in coef_[0] if i!=0], key=abs, reverse=True):
        samplefile = samplefiles[list(coef_[0]).index(coef)]
        dataID = samplefile.split('_')[0]
        if dataID in list(ann.keys()):
            annInfo =  ann[dataID]
        else:
            annInfo = ['NA','NA','NA']
        outf.write('{}\t{}\t{}\n'.format(dataID, coef, '\t'.join(annInfo)))

    outf.write('AUC = {}\n'.format(auc))
    outf.write('best_cvs = {}\n'.format(best_cvs))
    outf.write('selected_features = {}\n'.format(selected_features))

    return auc,selected_features

def main(args):
    '''
    Input arguments from command line.

    '''
    # read all parameters from arguments
    gxfile = args.expr # input gene symbols/refseqID
    name = args.name # output name
    exptype = args.exptype

    genome = args.genome # species
    symfile = args.sym
    annotation = args.annotation

    rp_hdf5 = args.histRP
    transform = args.transform

    maxsamples = args.maxsamples

    # TODO: symfile is the refseqID annotation file, change to gene symbol file?
    sym = gene_sym(symfile) # {"resfseq": {"gene_symbol", "chr"}}

    h5file = tables.open_file( rp_hdf5, driver="H5FD_CORE")
    samplenames  = getSampleNames_hdf5(h5file)
    z   = readregpotfiles(sym,genome,samplenames,h5file)
    h5file.close()
    
    if transform == 'log':
        z.x = logtransform(z.x)
    if transform == 'sqrt':
        z.x = sqrttransform(z.x)

    (genenames,z.y,train_index,test_index) = read_genelistOnly(sym, gxfile, z.index, exptype)

    sys.stdout.flush()
    print('Do regrssion with TARGET genes...')
    y = 1*( z.y == TARGET )

    x = z.x[:,:-5] # remove the last few columns: refseq, start, chr, etc...
    print("Adaptive lasso RP matrix shape...")
    print("{}\n".format(np.shape(x)))
    ann = dataset_annotation(annotation)
    try:
        auc,selected_features = adaptive_lasso(x,y,z.rpfiles,name,maxsamples,ann,genenames)
    except:
        sys.stderr.write("""\nERROR: bart2 exited with errors!
Please check whether you selected the correct species or uploaded the correct gene list!\n""")
        sys.exit(1)

    print("Adaptive lasso regression AUC score and selected features...")
    print(auc)
    print(selected_features)
    sys.stdout.flush()
 
if __name__ == "__main__":
    try:
        parser = argparse.ArgumentParser(description="""Regression of regulatory potential to gene expression changes.""")
        # related to input file type
        parser.add_argument( '-e','--expr', dest='expr', required = True, type = str, help = 'The related differential expression file')
        parser.add_argument( '--exptype', dest='exptype', required = True, choices=['Gene_Response','Gene_Only'], type = str, \
            help = 'Gene_Response includes 2 columns, one is the geneID, and the other is 1/0, 1 for target and 0 for un-target; \
                    Gene_Only includes 1 column, only the gene list of the targets. \
                    Only official gene symbol or refseqID are allowd for the geneID.')

        parser.add_argument( '-r','--historicalRP', dest='histRP', required = True, type = str, \
                help = 'The file with hdf5 format which contain the H3K27ac RP information')
        # transform method
        parser.add_argument('-t', '--transform', dest="transform", type=str, default='sqrt', choices=['sqrt', 'log'], required=True, \
                help='Use sqrt transform or log transform on RP ')
        parser.add_argument( '-a', dest='annotation', required=True,  type=str, help='The annotation file for each dataset' )

        parser.add_argument( '-m', dest='sym', type=str, required=True, help='refseqTSS is six columns: <chromosome name> <TSS> <TSS+1> <refseq:genesymbok> <score> <strand>')
        # parser.add_argument( '-m', dest='sym', type=str, required=True, help='genesymbolTSS is six columns: <chromosome name> <TSS-1> <TSS+1> <genesymbol> <score> <strand>')
        parser.add_argument( '-g','--genome', dest="genome", type=str, default='hg38', choices=['mm9','hg19','hg38','mm10'], required=False, help='genome')
        parser.add_argument( '--maxsamples', dest='maxsamples',  type=int, default=20, required=False, \
                help='Maximum number of samples to include in regression model.' )
        parser.add_argument( '-a', dest='annotation', required=True,  type=str, help='The annotation file for each dataset' )

        parser.add_argument( '-n','--name', dest='name',required = True, type = str, help = 'The prefix of the output names')

        args = parser.parse_args()
        main(args)

    except KeyboardInterrupt:
        sys.stderr.write("User interrupted me!\n")
        sys.exit(0)
