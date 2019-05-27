import os,sys
from types import SimpleNamespace

# import from package
from bart2 import OptValidator, ReadCount, RPRegress, EnhancerIdentifier, AUCcalc, StatTest, score_on_UDHS

script_dir = os.path.dirname(os.path.realpath(__file__))
ADAPTIVE_LASSO_MAXSAMPLES = 20 # TODO: should we fix it?

def bart(options):
    args = OptValidator.opt_validate(options)
 
    # create output directory
    try:
        os.makedirs(args.outdir, exist_ok=True)
    except:
        sys.stderr.write('Output directory: {} could not be created. \n'.format(args.outdir))
        sys.exit(1)
    sys.stdout.write("Output directory will be {} \n".format(args.outdir))
    sys.stdout.write("Output file prefix will be {} \n".format(args.ofilename))

    if args.species == 'hg38':
        sys.stdout.write("Start prediction on hg38...\n")
    elif args.species == 'mm10':
        sys.stdout.write("Start prediction on mm10...\n")
    sys.stdout.flush()

    # bart geneset [-h] <-i genelist.txt> [--refseq] <-s species> [-t target] [-p processes] [--outdir] [options]
    if args.subcommand_name == 'geneset':
        '''
        Use adaptive lasso regression to select H3K27ac samples.

        RPRegress parameters:
        species, rp matrix, gene file, refseq TSS, output directory, target gene method, adaptive lasso max sapmle numbers, log transform, square transform, gene symbol or refseqID, gene symbol or not, separate by chrom, sample description file
        '''
        sys.stdout.write("Do adaptive lasso to select informative H3K27ac samples...\n")
        sys.stdout.flush()

        sys.stdout.write("Generate parameters for regression step...\n")
        if options.refseq:
            rp_args = \
                SimpleNamespace(genome=args.species, \
                          histRP=args.rp, \
                          expr=args.infile, \
                          sym=args.tss, \
                          name=args.ofilename, \
                          maxsamples=ADAPTIVE_LASSO_MAXSAMPLES, \
                          transform='sqrt', \
                          exptype="Gene_Response", \
                          annotation=args.desc)
            RPRegress.main(rp_args)
        else:
            rp_args = \
                SimpleNamespace(genome=args.species, \
                          histRP=args.rp, \
                          expr=args.infile, \
                          sym=args.tss, \
                          name=args.ofilename, \
                          maxsamples=ADAPTIVE_LASSO_MAXSAMPLES, \
                          transform='sqrt', \
                          exptype="Gene_Only", \
                          annotation=args.desc)
            RPRegress.main(rp_args)

        regression_info = args.ofilename + '_adaptive_lasso_Info.txt'
        if not os.path.exists(regression_info):
            sys.stderr.write("Error: selecting samples from H3K27ac compendium! \n")
            sys.exit(1)

        '''
        Generate cis-regulatory profile based on adaptive lasso model weights multiply H3K27ac samples RPKM signal.

        EnhancerIdentifier parameters:
        selected samples file, gene file, output directory, species, UDHS, rpkm matrix
        '''
        sys.stdout.write("Generate cis-regulatory profile...\n")
        sys.stdout.flush()

        sys.stdout.write("Generate parameters for enhancer profile generation...\n")
        enhancer_args = \
                SimpleNamespace(samplefile=regression_info, \
                                name=args.ofilename, \
                                k27achdf5=args.rpkm)
        EnhancerIdentifier.main(enhancer_args)
        sys.stdout.flush()

        enhancer_profile = args.ofilename + '_enhancer_prediction_lasso.txt'
        if not os.path.exists(enhancer_profile):
            sys.stderr.write("Error: generating enhancer profile! \n")
            sys.exit(1)

        # get ranked score UDHS positions from enhancer profile
        positions = AUCcalc.get_position_list(enhancer_profile)

    # bart profile [-h] <-i ChIP-seq profile> <-f format> <-s species> [-t target] [-p processes] [--outdir] [options]
    elif args.subcommand_name == 'profile':
        sys.stdout.write('Start mapping the {} file...\n'.format(args.format.upper()))
        sys.stdout.flush()
        counting = ReadCount.read_count_on_DHS(args)
        # get ranked score UDHS positions from read count
        positions = sorted(counting.keys(),key=counting.get,reverse=True)

    elif args.subcommand_name == 'region':
        sys.stdout.write('Start mapping the bed score onto UDHS...\n')
        sys.stdout.flush()
        counting = score_on_UDHS.score_on_DHS(args)
        sys.stdout.write('Sorting scored UDHS...\n')
        sys.stdout.flush()
        positions = sorted(counting.keys(),key=counting.get,reverse=True)


    '''
    Start using revised BART on calculating the AUC score for each TF ChIP-seq dataset
    '''
    sys.stdout.write('BART Prediction starts...\n\nRank all DHS...\n')
    sys.stdout.flush()

    tf_aucs, tf_index = AUCcalc.cal_auc(args, positions)
    sys.stdout.flush()

    stat_file = args.ofilename + '_bart_results.txt'
    StatTest.stat_test(tf_aucs, tf_index, stat_file, args.normfile)
    sys.stdout.flush()
    sys.stdout.write("Congratulations! BART job finished successfully!\n")

