#!/usr/bin/env python3
#-*- coding: utf-8 -*-

import os, re, sys, argparse
import numpy as np
import pandas as pd
import subprocess
import configparser
import warnings
warnings.filterwarnings('ignore')
import logging
from gc_correct import GC_correct

bindir = os.path.abspath(os.path.dirname(__file__))
pat1 = re.compile(r'^\s+$')

def keep_same_order(infile, bed, outfile):
    sample_info = dict()
    with open(infile) as IN:
        for line in IN:
            if line.startswith('chr\tstart'):
                header = line.strip().split('\t')
                continue
            else:
                items = dict(zip(header, line.strip().split()))
                key = f'{items["chr"]}\t{items["start"]}\t{items["end"]}\t{items["gene"]}\t{items["gc"]}\t{items["strand"]}\t{items["exon"]}'
                value = f'{items["mean_depth_target"]}\t{items["mean_depth_autochr"]}\t{items["normalized_depth"]}\t{items["normalized_depth_gc"]}'
                sample_info[key] = value

    with open(bed) as BED, open(outfile, 'w') as OUT:
        OUT.write('\t'.join(header) + '\n')
        for line in BED:
            chr, start, end, gene, gc, strand, transcript, exon = line.strip().split('\t')
            key = f'{chr}\t{start}\t{end}\t{gene}\t{gc}\t{strand}\t{exon}'
            if key in sample_info:
                OUT.write(line.strip() + '\t' + sample_info[key] + '\n')
            else:
                sys.exit(f'{key} not found in {bed}')

def cal_mean_depth(bam, outfile, bedfile, config, gc_correct=False, rm_tmp=False):
    outdir = os.path.dirname(outfile)
    sample = os.path.basename(outfile).split('.')[0]
    samtools = config.get('software', 'samtools')
    bedtools = config.get('software', 'bedtools')
    warning_depth_threshold = int(config.get('parameter', 'warning_depth_threshold'))
    fail_depth_threshold = int(config.get('parameter', 'fail_depth_threshold'))

    bedcov_cmd = f'{samtools} bedcov {bedfile} {bam} > {outfile}.bedcov'
    subprocess.call(bedcov_cmd, shell=True)
    bedcov_df = pd.read_csv(f'{outfile}.bedcov', sep='\t', names=['chr', 'start', 'end', 'gene', 'gc', 'strand', 'transcript', 'exon', 'total_depth'])
    bedcov_auto = bedcov_df.loc[((bedcov_df['chr']!='chrX')&(bedcov_df['chr']!='chrY')),:] # skip sex chr
    bedcov_auto['mean_depth_target'] = bedcov_auto.apply(lambda x: int(x['total_depth']/(x['end'] - x['start'])), axis=1) # cal mean depth
    bedcov_auto['mean_depth_autochr'] = int(bedcov_auto['mean_depth_target'].mean()) # add mean depth
    bedcov_auto['normalized_depth'] = np.round(bedcov_auto['mean_depth_target']/bedcov_auto['mean_depth_autochr'], 2) # cal normalized depth
    bedcov_auto.loc[:, 'gc'] = bedcov_auto['gc'].apply(lambda x: format(x, '.2f'))
    if gc_correct:
        bedcov_auto = GC_correct(bedcov_auto).rolling_median()
        bedcov_auto[['chr', 'start', 'end', 'gene', 'gc', 'strand', 'transcript', 'exon', 'mean_depth_target', 'mean_depth_autochr', 'normalized_depth', 'normalized_depth_gc']].to_csv(outfile + '.tmp', sep='\t', index=False)
    else:
        bedcov_auto[['chr', 'start', 'end', 'gene', 'gc', 'strand', 'transcript', 'exon', 'mean_depth_target', 'mean_depth_autochr', 'normalized_depth']].to_csv(outfile + '.tmp', sep='\t', index=False)
    keep_same_order(outfile + '.tmp', bedfile, outfile)


    # # warnings
    # bedcov_auto = bedcov_auto[~((bedcov_auto['gene']=='PDGFRA')&(bedcov_auto['start']==55106199))] # no probe target
    # if (bedcov_auto['mean_depth_target'] < fail_depth_threshold).any():
    #     logging.error('Exist target that its depth is ultra low: <50x, please check fail.txt\n')
    #     bedcov_auto[bedcov_auto['mean_depth_target'] < fail_depth_threshold].to_csv(f'{outdir}/{sample}.fail.txt', sep='\t', index=False)
    # elif (bedcov_auto['mean_depth_target'] < warning_depth_threshold).any():
    #     logging.warning('Exist target that its depth is low: 50-100x, please check warning.txt\n')
    #     bedcov_auto[bedcov_auto['mean_depth_target'] < warning_depth_threshold].to_csv(f'{outdir}/{sample}.warning.txt', sep='\t', index=False)

    if rm_tmp:
        subprocess.call(f'rm {outfile}.bedcov', shell=True)
        subprocess.call(f'rm {outfile}.tmp', shell=True)


def main():
    parser = argparse.ArgumentParser(usage='\nInfer Germline CNV (LGR) Using WBC-PoN Method', description=__doc__)
    parser.add_argument('-i', '--inbam', help='bam file', dest='inbam', required=True)
    parser.add_argument('-s', '--sample', help='sample name', dest='sample', required=True)
    parser.add_argument('-o', '--outdir', help='outdir', dest='outdir', required=True)
    parser.add_argument('-b', '--bed', help='bed file', dest='bed', required=True)
    parser.add_argument('--gc_correct', help='Do gc correction', dest='gc_correct', action='store_true', default=False)
    parser.add_argument('-c', '--config', help='config file', dest='config')
    parser.add_argument('--samtools', help='samtools dir', dest='samtools', default='/share/Onc_Soft_DB/software/samtools-1.4/samtools')
    args = parser.parse_args()

    configfile = args.config if args.config else f'{bindir}/../config/config.ini'
    config = configparser.ConfigParser()
    config.read(configfile, encoding='utf-8')

    # cal mean depth of each target
    cal_mean_depth(args.inbam, f'{args.outdir}/{args.sample}.normalized_depth.txt', args.bed, config, args.gc_correct, rm_tmp=True)

if __name__ == "__main__":
    main()
