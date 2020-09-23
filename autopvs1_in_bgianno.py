# -*- coding:utf-8 -*-
from multiprocessing import Pool, cpu_count
from autopvs1 import AutoPVS1
from functools import reduce
import pandas as pd
import optparse
import pyfaidx
import yaml
import sys
import os


def print_usage(option, opt, value, parser):
    usage_message = r"""
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
    python3 autopvs1_in_bgianno.py -i test.bed -o test.autopvs1.tsv -p 1
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
    """
    print(usage_message)
    sys.exit()


def yaml_read(yaml_file):
    with open(yaml_file, 'r') as y:
        yaml_dic = yaml.load(y, Loader=yaml.FullLoader)
    return yaml_dic


def bgi_anno_2_vcf_format(df, fa):
    df.reset_index(drop=True, inplace=True)
    df['#Chr'] = df['#Chr'].astype('str')
    if len(df[df['#Chr'].str.startswith('chr')]):
        df['#CHROM'] = df['#Chr']
    else:
        df['#CHROM'] = 'chr' + df['#Chr']
    df.loc[df['#CHROM'] == 'chrMT', '#CHROM'] = 'chrM_NC_012920.1'
    df['ID'] = '.'
    df['QUAL'] = '.'
    df['FILTER'] = '.'
    df['INFO'] = '.'
    df['MuType'] = 'delins'
    df.loc[df['Ref'] == '.', 'MuType'] = 'ins'
    df.loc[df['Call'] == '.', 'MuType'] = 'del'
    df.loc[(df['Ref'].map(len) == 1) & (df['Call'].map(len) == 1) & (df['Ref'] != '.') & (df['Call'] != '.'), 'MuType'] = 'snp'
    df['POS'] = df['Stop']
    df.loc[df['MuType'] == 'del', 'POS'] = df.loc[df['MuType'] == 'del', 'Start']
    df.loc[df['MuType'] == 'delins', 'POS'] = df.loc[df['MuType'] == 'delins', 'Start']
    df['REF'] = df['Ref']
    df['ALT'] = df['Call']
    for i in range(df.shape[0]):
        if df.loc[i, 'MuType'] == 'ins':
            base = str(fa.get_seq(df.loc[i, '#CHROM'], df.loc[i, 'POS'], df.loc[i, 'POS'])).upper()
            df.loc[i, 'REF'] = base
            df.loc[i, 'ALT'] = base + df.loc[i, 'ALT']
        elif df.loc[i, 'MuType'] == 'del':
            base = str(fa.get_seq(df.loc[i, '#CHROM'], df.loc[i, 'POS'], df.loc[i, 'POS'])).upper()
            df.loc[i, 'ALT'] = base
            df.loc[i, 'REF'] = base + df.loc[i, 'REF']
        elif df.loc[i, 'MuType'] == 'delins':
            base = str(fa.get_seq(df.loc[i, '#CHROM'], df.loc[i, 'POS'], df.loc[i, 'POS'])).upper()
            df.loc[i, 'REF'] = base + df.loc[i, 'REF']
            df.loc[i, 'ALT'] = base + df.loc[i, 'ALT']
        else:
            pass
    a = df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']].copy()
    a.sort_values(by=['#CHROM', 'POS'], ascending=True, inplace=True)
    b = df[['#Chr', 'Start', 'Stop', 'Ref', 'Call', '#CHROM', 'POS', 'REF', 'ALT']].copy()
    df.drop(columns=['ID', 'QUAL', 'FILTER', 'INFO', 'MuType'], inplace=True)
    return a, b, df


def vcf_format_2_bgi_anno(in_df):
    df = in_df.copy()
    df['MuType'] = '.'
    df.loc[(df['REF'].str.len() == 1) & (df['ALT'].str.len() == 1), 'MuType'] = 'snv'
    df.loc[(df['REF'].str.len() == 1) & (df['ALT'].str.len() > 1), 'MuType'] = 'ins'
    df.loc[(df['REF'].str.len() > 1) & (df['ALT'].str.len() == 1), 'MuType'] = 'del'
    df.loc[(df['REF'].str.len() > 1) & (df['ALT'].str.len() > 1) & (df['REF'].str[0] == df['ALT'].str[0]), 'MuType'] = 'delins_eq'
    df.loc[(df['REF'].str.len() > 1) & (df['ALT'].str.len() > 1) & (df['REF'].str[0] != df['ALT'].str[0]), 'MuType'] = 'delins_neq'
    df['#CHROM'] = df['#CHROM'].astype('str')
    if len(df[df['#CHROM'].str.startswith('chr')]):
        df['#Chr'] = df['#CHROM']
    else:
        df['#Chr'] = 'chr' + df['#CHROM']
    df.loc[df['#CHROM'] == 'chrM_NC_012920.1', '#Chr'] = 'chrMT'
    df['Start'] = df['POS'].astype(int) - 1
    df['Stop'] = df['POS'].astype(int)
    df['Ref'] = df['REF']
    df['Call'] = df['ALT']
    df.loc[df['MuType'] == 'ins', 'Ref'] = '.'
    df.loc[df['MuType'] == 'ins', 'Call'] = df.loc[df['MuType'] == 'ins', 'ALT'].str[1:]
    df.loc[df['MuType'] == 'del', 'Call'] = '.'
    df.loc[df['MuType'] == 'del', 'Ref'] = df.loc[df['MuType'] == 'del', 'REF'].str[1:]
    df.loc[df['MuType'] == 'ins', 'Start'] = df.loc[df['MuType'] == 'ins', 'Stop']
    df.loc[df['MuType'] == 'del', 'Start'] = df.loc[df['MuType'] == 'del', 'Stop']
    df.loc[df['MuType'] == 'del', 'Stop'] = df.loc[df['MuType'] == 'del', 'Stop'] + df.loc[df['MuType'] == 'del', 'Ref'].str.len()
    df.loc[df['MuType'] == 'delins_eq', 'Stop'] = df.loc[df['MuType'] == 'delins_eq', 'Start'] + df.loc[df['MuType'] == 'delins_eq', 'Ref'].str.len()
    df.loc[df['MuType'] == 'delins_eq', 'Start'] = df.loc[df['MuType'] == 'delins_eq', 'POS']
    df.loc[df['MuType'] == 'delins_eq', 'Ref'] = df.loc[df['MuType'] == 'delins_eq', 'REF'].str[1:]
    df.loc[df['MuType'] == 'delins_eq', 'Call'] = df.loc[df['MuType'] == 'delins_eq', 'ALT'].str[1:]
    df.loc[df['MuType'] == 'delins_neq', 'Stop'] = df.loc[df['MuType'] == 'delins_neq', 'Start'] + df.loc[df['MuType'] == 'delins_neq', 'Ref'].str.len()
    a = df[['#CHROM', 'POS', 'REF', 'ALT', '#Chr', 'Start', 'Stop', 'Ref', 'Call']].copy()
    df.drop(columns=['MuType'], inplace=True)
    return a, df


def autopvs1_anno(mut, trans):
    try:
        autopvs1_res = AutoPVS1(mut, trans)
        hgvs_c, strength, strength_adj, criterion = str(autopvs1_res.hgvs_c), str(autopvs1_res.pvs1.strength_raw), str(autopvs1_res.pvs1.strength), str(autopvs1_res.pvs1.criterion)
    except:
        hgvs_c, strength, strength_adj, criterion = '.', '.', '.', '.'
    return ';'.join([hgvs_c, strength, strength_adj, criterion])


def run_autopvs1_anno(df):
    df.reset_index(drop=True, inplace=True)
    df['AutoPVS1 input'] = df['#CHROM'] + '-' + df['POS'].astype('str') + '-' + df['REF'] + '-' + df['ALT']
    df['AutoPVS1 result'] = df.apply(lambda x: autopvs1_anno(x['AutoPVS1 input'], x['Transcript']), axis=1)
    aupvs1_res = df['AutoPVS1 result'].str.split(';', expand=True)
    aupvs1_res.columns = ['AutoPVS1 Hgvs_c', 'AutoPVS1 Strength', 'AutoPVS1 Adjusted Strength', 'AutoPVS1 Criterion']
    df = df.join(aupvs1_res)
    df.replace({'AutoPVS1 Strength': 'Strength.', 'AutoPVS1 Adjusted Strength': 'Strength.'}, {'AutoPVS1 Strength': '', 'AutoPVS1 Adjusted Strength': ''}, regex=True, inplace=True)
    df.drop(columns=['#CHROM', 'POS', 'REF', 'ALT', 'AutoPVS1 input', 'AutoPVS1 result'], inplace=True)
    return df


def split_df(df, split_num):
    df_list = list()
    step = round(df.shape[0]/split_num)
    for i in range(split_num):
        if i == 0:
            df_list.append(df.loc[0: step-1])
        elif i == split_num-1:
            df_list.append(df.loc[step*i:])
        else:
            df_list.append(df.loc[step*i:step*(i+1)-1])
    return df_list


path = os.path.split(os.path.realpath(__file__))[0]
if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('-u', '--usage', help='print more info on how to use this script', action="callback", callback=print_usage)
    parser.add_option('-i', '--in', dest='in_file', help='input file', default=None, metavar='file')
    parser.add_option('--skip_rows', dest='skip_rows', default=0, type=int)
    parser.add_option('-o', '--out', dest='out', help='output file', default=None, metavar='string')
    parser.add_option('-p', '--process', dest='process', help='process num', default=1, type=int)
    parser.add_option('--pwd', dest='pwd', default=path, metavar='string')
    parser.add_option('-c', '--config', dest='config', default=os.path.join(path, 'etc', 'autopvs1.yaml'), metavar='file')
    (opts, args) = parser.parse_args()
    in_file = opts.in_file
    skip_rows = opts.skip_rows
    out = opts.out
    pwd = opts.pwd
    config = opts.config
    process = opts.process
    process_num = min(process, cpu_count())
    config_dic = yaml_read(config)
    reference = os.path.join(pwd, config_dic['reference'])
    lof_genes_file = os.path.join(pwd, config_dic['lof_genes'])
    if not os.path.exists(config):
        print(config + ' is not exist')
        sys.exit(1)
    if not os.path.exists(reference):
        print(reference + ' is not exist')
        sys.exit(1)
    if not os.path.exists(lof_genes_file):
        print(lof_genes_file + ' is not exist')
        sys.exit(1)
    fa_read = pyfaidx.Fasta(reference)
    lof_genes = pd.read_csv(lof_genes_file, sep='\t', header=None)
    null_list = ['splice-3', 'splice-5', 'init-loss', 'alt-start', 'frameshift', 'nonsense', 'stop-gain', 'span']
    anno_df = pd.read_csv(in_file, sep='\t', dtype={'#Chr': str}, low_memory=False, skiprows=range(skip_rows))
    anno_df_not_lof = anno_df[~((anno_df['Gene Symbol'].isin(lof_genes[0].values)) & (anno_df['Function'].isin(null_list)) & (anno_df['Ref'] != anno_df['Call']))].copy()
    if not anno_df_not_lof.empty:
        anno_df_not_lof['AutoPVS1 Hgvs_c'] = '.'
        anno_df_not_lof['AutoPVS1 Strength'] = '.'
        anno_df_not_lof['AutoPVS1 Adjusted Strength'] = '.'
        anno_df_not_lof['AutoPVS1 Criterion'] = '.'
    else:
        anno_df_not_lof['AutoPVS1 Hgvs_c'] = []
        anno_df_not_lof['AutoPVS1 Strength'] = []
        anno_df_not_lof['AutoPVS1 Adjusted Strength'] = []
        anno_df_not_lof['AutoPVS1 Criterion'] = []
    anno_df_lof = anno_df[(anno_df['Gene Symbol'].isin(lof_genes[0].values)) & (anno_df['Function'].isin(null_list)) & (anno_df['Ref'] != anno_df['Call'])].copy()
    if not anno_df_lof.empty:
        _1, _2, anno_df_lof = bgi_anno_2_vcf_format(anno_df_lof, fa_read)
        anno_df_lof_list = split_df(anno_df_lof, process_num)
        with Pool(process_num) as pool:
            df_pvs1_list = pool.map(run_autopvs1_anno, anno_df_lof_list)
        df_pvs1_merge = reduce(lambda x, y: x.append(y), df_pvs1_list)
        df_pvs1_merge = df_pvs1_merge.append(anno_df_not_lof, sort=False)
        df_pvs1_merge.to_csv(out, sep='\t', index=False)
    else:
        anno_df_not_lof.to_csv(out, sep='\t', index=False)
