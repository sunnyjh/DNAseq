#!/usr/bin/env python3
#-*- coding: utf-8 -*-

__author__ = 'Zhang Yikai'
__mail__ = 'zhangyk4498@berryoncology.com'
__usage__= 'create ref matrix by PoN file'

import os, re, sys, argparse
import glob
import logging

bindir = os.path.abspath(os.path.dirname(__file__))
part1 = re.compile(r'^\s+$')

class make_ref_matrix(object):
    def __init__(self,  refdir, outfile) -> None:
        self.refdir = refdir
        self.outfile = outfile

    def get_target(self, cnn_file):
        target_dict = {}
        with open(cnn_file, 'r') as IN:
            for line in IN:
                if line.startswith('chromosome'):
                    continue
                line = line.strip('\n').split('\t')
                key = f'{line[0]}.{line[1]}.{line[2]}.{line[3]}'
                target_dict[key] = []
        return target_dict

    def get_samples(self, cnn_file, samples, target_dict):
        sample = os.path.basename(cnn_file).split('.')[0]
        samples.append(sample)
        with open(cnn_file, 'r') as IN:
            for line in IN:
                if line.startswith('chromosome'):
                    continue
                chrom, start, end, gene, log2, depth, gc, spread = line.strip('\n').split('\t')
                key = f'{chrom}.{start}.{end}.{gene}'
                values = f'{log2};{depth}'
                if key in target_dict.keys():
                    target_dict[key].append(values)
        return samples, target_dict

    def Running(self):
        '''
        get ref matrix
        '''

        cnn_list = sorted(glob.glob(f'{self.refdir}/*/*.cnr'))
        cnn_list_len = len(cnn_list)
        print(f'Find {cnn_list_len} refs\n')

        target_dict = self.get_target(cnn_list[0])
        samples = []
        for file in cnn_list:
            samples, target_dict = self.get_samples(file, samples, target_dict)
        
        samples_list = "\t".join(samples)
        with open(self.outfile, 'w') as OUT:
            OUT.write(f'Target\t{samples_list}\n')
            for key in target_dict.keys():
                values = "\t".join(target_dict[key])
                OUT.write(f'{key}\t{values}\n')

def main():
    parser = argparse.ArgumentParser(usage=f'{__usage__}',
                description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter,
                epilog='author:\t{0}\nmail:\t{1}'.format(__author__, __mail__))
    parser.add_argument('-refdir', help='the ref PoN directory...', dest='refdir', required=True)
    parser.add_argument('-outfile',  help='output file...', dest='outfile', required=True)
    args = parser.parse_args()

    # set the logging
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s")
    logging.info('step1: running ,make ref matrix')
    run = make_ref_matrix(args.refdir, args.outfile)
    run.Running()

if __name__ == "__main__":
    main()
