import os
import sys
import time
import pandas as pd


def write_auc(args, tf_auc, tf_dict):
    outf = open(args.ofilename + '_auc.txt', 'w')
    for key in sorted(tf_auc.keys(), key=tf_auc.get, reverse=True):
        tf_file_name = tf_dict[str(key)]
        outf.write('{}\tAUC = {:.3f}\n'.format(
            tf_file_name, tf_auc[key]))
    outf.close()
    sys.stdout.write(
        '\n--ROC-AUC calculation finished!\n--Results saved in file: {}\n'.format(args.ofilename + '_auc.txt'))


def write_adaptive_lasso(args, main_df, additional_lines):
    outf = open('{}_adaptive_lasso_Info.txt'.format(args.ofilename), 'w')

    for i in range(0, len(main_df.index)):
        out_line = main_df.loc[i, ].values.flatten().tolist()
        out_line = list(map(str, out_line))
        outf.write("\t".join(out_line) + '\n')

    outf.write('AUC = {}\n'.format(additional_lines[0]))
    outf.write('best_cvs = {}\n'.format(additional_lines[1]))
    outf.write('selected_features = {}\n'.format(additional_lines[2]))

    outf.close()

def write_udhs_prediction(args, df):
    outf = open(args.ofilename+'_enhancer_prediction_lasso.txt','w')
    outf.write("\t".join(list(df.columns)) + '\n')
    for i in range(0, len(df.index)):
        out_line = df.loc[i, ].values.flatten().tolist()
        out_line[4] = float(out_line[4])
        outf.write('{}\t{}\t{}\t{}\t{:.2f}\n'.format(out_line[0], out_line[1], out_line[2], out_line[3], out_line[4]))

    outf.close()

def write_bart_result(args, df):
    outf = open(args.ofilename + '_bart_results.txt', 'w')
    outf.write('TR\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('statistic','pvalue','zscore','max_auc','re_rank','irwin_hall_pvalue'))
    for i in df.index:
        outf.write('{}\t{:.3f}\t{:.3e}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3e}\n'.format(i,df.loc[i]['score'],df.loc[i]['pvalue'],df.loc[i]['zscore'],df.loc[i]['max_auc'],df.loc[i]['rank_avg_z_p_a'],df.loc[i]['rank_avg_z_p_a_irwinhall_pvalue']))
    sys.stdout.write('--Standardization finished!\n--Ranked TRs saved in file: {}\n'.format(args.ofilename + '_bart_results.txt'))

