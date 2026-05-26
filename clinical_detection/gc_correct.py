#!/usr/bin/env python3
#-*- coding: utf-8 -*-

import os, re, sys, argparse
import configparser
import numpy as np
import pandas as pd
import math
import subprocess
import warnings
warnings.filterwarnings('ignore')

bindir = os.path.abspath(os.path.dirname(__file__))
pat1 = re.compile(r'^\s+$')

bedtools = '/share/Onc_Soft_DB/software/bedtools2/2.26.0/bin/bedtools'

class GC_correct():
    def __init__(self, in_df):
        self.in_df = in_df
        self.frac = 0.1
        self.min_wing = 3

    def rolling_median(self):
        gc = self.in_df['gc']
        order = np.argsort(gc, kind='mergesort')
        self.in_df = self.in_df.iloc[order]
        depth = self.in_df['normalized_depth']
        wing = int(math.ceil(len(depth) * self.frac * 0.5))
        wing = max(wing, self.min_wing)
        wing = min(wing, len(depth) - 1)
        signal = np.concatenate((depth[wing-1::-1], depth, depth[:-wing-1:-1]))
        signal = pd.Series(signal)
        rolled = signal.rolling(2 * wing + 1, 1, center=True).median()
        biases = np.asfarray(rolled[wing:-wing])
        self.in_df['normalized_depth_gc'] = round(self.in_df['normalized_depth'] / biases, 2)
        return self.in_df

    def lowess(self):
        pass


def main():
    parser = argparse.ArgumentParser(usage='\nInfer Germline CNV (LGR) Using WBC-PoN Method', description=__doc__)
    parser.add_argument('-i', '--input', help='input depth file', dest='input', required=True)
    parser.add_argument('-o', '--output', help='output file', dest='output', required=True)
    parser.add_argument('-b', '--bed', help='bed file', dest='bed')
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep='\t')
    if 'gc' not in df.columns:
        subprocess.call(f'{bedtools} intersect -a {args.bed} -b {args.input} -wa -wb | cut -f1-8,13-15|sed "1i chr\tstart\tend\tgene\tgc\tstrand\ttranscript\texon\tmean_depth_target\tmean_depth_autochr\tnormalized_depth" > {args.output}.tmp', shell=True)
        df = pd.read_csv(f'{args.output}.tmp', sep='\t')
    df.loc[:, 'gc'] = df['gc'].apply(lambda x: format(x, '.2f'))
    do_gc = GC_correct(df)
    df_gc = do_gc.rolling_median()
    df_gc.to_csv(args.output, sep='\t', index=False)
    subprocess.call(f'sort -k1,1 -V {args.output} > {args.output}.tmp && mv {args.output}.tmp {args.output}', shell=True)


if __name__ == "__main__":
    main()
