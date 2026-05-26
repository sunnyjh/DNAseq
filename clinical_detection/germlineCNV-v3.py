#!/usr/bin/env python3
#-*- coding: utf-8 -*-

import os, re, sys, argparse
import configparser
import glob
import warnings
warnings.filterwarnings('ignore')
from cal_depth import cal_mean_depth
from deal_PoN import PoN
from infer_CNV import CNV
from merge_CNV import mergeCNV
from get_Report import getReport, getExon_num
import monitor_ultra_low_depth
import logging
import subprocess

bindir = os.path.abspath(os.path.dirname(__file__))
pat1 = re.compile(r'^\s+$')

def collect_PoN_bak(PoN_dir=None, panel=None, seqtype=None, old=False, suffix='normalized_depth.txt'):
    if PoN_dir:
        assert os.path.exists(PoN_dir)
        return glob.glob(f'{PoN_dir}/*{suffix}')
    elif panel and seqtype and not old:
        if seqtype != "T7" and seqtype != "t7":
            seqtype = 'NOVA'
        else:
            seqtype = 'T7'
        if panel == '2':
            panel = '31'
        return glob.glob(f'{bindir}/../database/PoN/{panel}/{seqtype}/*{suffix}')
    elif panel and old:
        if seqtype != "T7" and seqtype != "t7":
            seqtype = 'NOVA'
        if panel == '2':
            panel = '31'
        return glob.glob(f'{bindir}/../database/PoN/{panel}/old/*{suffix}')

def collect_PoN(PoN_dir=None, panel=None, seqtype=None, lib_version=None, cons=False, suffix='normalized_depth.txt'):
    if PoN_dir:
        if len(glob.glob(f'{PoN_dir}/*{suffix}')) > 0:
            return glob.glob(f'{PoN_dir}/*{suffix}')
        else:
            if lib_version == 'vtp.6.1.4.6:vu.3.2.4.3':
                if seqtype != "T7" and seqtype != "t7":
                    seqtype = 'NOVA'
                else:
                    seqtype = 'T7'
                if panel == '2':
                    panel = '31'
                return glob.glob(f'{PoN_dir}/cons_bam/{panel}/{seqtype}/*{suffix}')
            else:
                if seqtype != "T7" and seqtype != "t7":
                    seqtype = 'NOVA'
                if panel == '2':
                    panel = '31'
                return glob.glob(f'{PoN_dir}/cons_bam/{panel}/old/*{suffix}')
    else:
        if not cons:
            if lib_version == 'vtp.6.1.4.6:vu.3.2.4.3':
                if seqtype != "T7" and seqtype != "t7":
                    seqtype = 'NOVA'
                else:
                    seqtype = 'T7'
                if panel == '2':
                    panel = '31'
                return glob.glob(f'{bindir}/../database/PoN/sort_bam/{panel}/{seqtype}/*{suffix}')
            else:
                if seqtype != "T7" and seqtype != "t7":
                    seqtype = 'NOVA'
                if panel == '2':
                    panel = '31'
                return glob.glob(f'{bindir}/../database/PoN/sort_bam/{panel}/old/*{suffix}')

        else:
            if lib_version == 'vtp.6.1.4.6:vu.3.2.4.3':
                if seqtype != "T7" and seqtype != "t7":
                    seqtype = 'NOVA'
                else:
                    seqtype = 'T7'
                if panel == '2':
                    panel = '31'
                return glob.glob(f'{bindir}/../database/PoN/cons_bam/{panel}/{seqtype}/*{suffix}')
            else:
                if seqtype != "T7" and seqtype != "t7":
                    seqtype = 'NOVA'
                if panel == '2':
                    panel = '31'
                return glob.glob(f'{bindir}/../database/PoN/cons_bam/{panel}/old/*{suffix}')

def parse_config(configfile):
    config = configparser.ConfigParser()
    config.read(configfile, encoding='utf-8')
    return config

def get_PlatForm(infile):
    seqtype_dict = {}
    lib_version_dict = {}
    with open(infile, 'r') as IN:
        for line in IN:
            sample, ol_check_process, seqtype = line.strip('\n').split('\t')
            if sample not in seqtype_dict:
                seqtype_dict[sample] = seqtype
            if sample not in lib_version_dict:
                lib_version_dict[sample] = ol_check_process
    return seqtype_dict, lib_version_dict

def main():
    parser = argparse.ArgumentParser(usage='\nInfer Germline CNV (LGR) Using WBC-PoN Method', description=__doc__)
    parser.add_argument('-i', '--inbam', help='bam file', dest='inbam', required=True)
    parser.add_argument('-s', '--sample', help='sample name', dest='sample', required=True)
    parser.add_argument('-o', '--outdir', help='outdir', dest='outdir', required=True)
    parser.add_argument('-p', '--panel', help='panel:[2, 31, 102, wes]', dest='panel', required=True)
    parser.add_argument('-m', '--method', help='matchscore method:[delta, corr]', dest='method', default='corr', choices=['corr', 'delta'])
    parser.add_argument('-x', '--select_method', help='select method:[mix, fix]', dest='select_method', default='mix', choices=['mix', 'fix'])
    parser.add_argument('-g', '--gc_correct', help='Do gc correction', dest='gc_correct', action='store_true', default=False)
    parser.add_argument('-b', '--bed', help='bed file', dest='bed')
    parser.add_argument('-r', '--refGeneWithVer', help='refGeneWithVer file', dest='refGeneWithVer')
    parser.add_argument('-d', '--PoN_dir', help='PoN dir', dest='PoN_dir', default='/share/Onc_Soft_DB/database/germline/LGR/PoN')
    parser.add_argument('-c', '--config', help='config file', dest='config', default=os.path.abspath(f'{bindir}/../config/config.ini'))
    parser.add_argument('--cons', help='Use cons bam PoN dir', dest='cons', action="store_true", default=False)
    parser.add_argument('--seqtype', help='Sample sequence type', dest='seqtype')
    args = parser.parse_args()

    # set the logging
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s")
    
    # parse config
    if args.panel == "wes":
        args.config = f'{bindir}/../config/config_wes.ini'
    config = parse_config(args.config)

    args.panel = '31' if args.panel == '2' else args.panel

    # get bed file
    args.bed = args.bed if args.bed else f'{bindir}/../database/Anno/{args.panel}gene.anno.bed'
	
    # get refGeneWithVer file
    args.refGeneWithVer = args.refGeneWithVer if args.refGeneWithVer else config.get('database', 'hg19_refGeneWithVer')

    # check input
    assert os.path.exists(args.inbam)
    assert os.path.exists(args.config)
    assert os.path.exists(args.bed)
    assert os.path.exists(args.refGeneWithVer)
    os.makedirs(args.outdir, exist_ok=True)

    # files
    sample_depth_file = f'{args.outdir}/{args.sample}.normalized_depth.txt'
    refmatrix = f'{args.outdir}/ref_matrix.xls'
    # ratio_file = f'{args.outdir}/{args.sample}.ratio.xls'
    anno_file = f'{args.outdir}/{args.sample}.anno.xls'
    anno_file_filter_UTR = f'{args.outdir}/{args.sample}.anno.filter.UTR.xls'
    anno_file_cutoff = f'{args.outdir}/{args.sample}.anno.cutoff.xls'
    exon_level_cnv = f'{args.outdir}/{args.sample}.exon.cnv.xls'
    # gene_level_cnv = f'{args.outdir}/{args.sample}.gene.cnv.xls'
    report_file = f'{args.outdir}/{args.sample}.rpt.cnv.xls'
    fig_cn_pdf = f'{args.outdir}/{args.sample}.CNV.pdf'
    fig_cn_png = f'{args.outdir}/{args.sample}.CNV.png'
    
    # gc_correct
    gc_correct = '--gc_correct' if args.gc_correct else ''

    # collect PoN samples
    if args.seqtype:
        seqtype_dict, lib_version_dict = get_PlatForm(args.seqtype)
        seqtypeinfo=seqtype_dict[args.sample]
        lib_version = lib_version_dict[args.sample]
    else:
        try:
            seqtype_dict, lib_version_dict = get_PlatForm(f'{args.outdir}/../../../../ExpVersion_PlatForm')
        except:
            seqtype_dict, lib_version_dict = get_PlatForm(f'{args.outdir}/../../../../../ExpVersion_PlatForm')
            # seqtype_dict, lib_version_dict = get_PlatForm(f'{args.outdir}/ExpVersion_PlatForm')
        seqtypeinfo=seqtype_dict[args.sample]
        lib_version = lib_version_dict[args.sample]
    PoN_list = collect_PoN(PoN_dir=args.PoN_dir, panel=args.panel, seqtype=seqtypeinfo, lib_version=lib_version, cons=args.cons) if args.PoN_dir else collect_PoN(panel=args.panel, seqtype=seqtypeinfo, lib_version=lib_version, cons=args.cons)

    # step1: cal mean depth of each target
    if args.inbam.endswith('.bam'):
        logging.info('step1: cal mean depth of each target')
        print(f'python3 {bindir}/cal_depth.py -i {args.inbam} -s {args.sample} -o {args.outdir} -b {args.bed} {gc_correct}\n')
        cal_mean_depth(args.inbam, sample_depth_file, args.bed, config, gc_correct=args.gc_correct, rm_tmp=True)
        monitor_ultra_low_depth.monitor_low_depth(sample_depth_file, f'{args.outdir}/warnings.txt', config)
    elif args.inbam.endswith('.normalized_depth.txt'):
        if args.gc_correct:
            logging.info('step1: do gc correct for normalized depth file')
            print(f'python3 {bindir}/gc_correct.py -i {args.inbam} -o {sample_depth_file}.tmp -b {args.bed} && mv {sample_depth_file}.tmp {sample_depth_file}\n')
            subprocess.call(f'python3 {bindir}/gc_correct.py -i {args.inbam} -o {sample_depth_file}.tmp -b {args.bed} && mv {sample_depth_file}.tmp {sample_depth_file}', shell=True)
            monitor_ultra_low_depth.monitor_low_depth(sample_depth_file, f'{args.outdir}/warnings.txt', config)
        else:
            logging.info('step1: intersect normalized depth file to outdir')
            print(f'{config.get("software", "bedtools")} intersect -a {args.bed} -b {args.inbam} -wa -wb | cut -f1-8,13-15|sed "1i chr\\tstart\\tend\\tgene\\tgc\\tstrand\\ttranscript\\texon\\tmean_depth_target\\tmean_depth_autochr\\tnormalized_depth" > {sample_depth_file}.tmp && mv {sample_depth_file}.tmp {sample_depth_file}\n')
            subprocess.call(f'{config.get("software", "bedtools")} intersect -a {args.bed} -b {args.inbam} -wa -wb | cut -f1-8,13-15|sed "1i chr\tstart\tend\tgene\tgc\tstrand\ttranscript\texon\tmean_depth_target\tmean_depth_autochr\tnormalized_depth" > {sample_depth_file}.tmp && mv {sample_depth_file}.tmp {sample_depth_file}', shell=True)
            monitor_ultra_low_depth.monitor_low_depth(sample_depth_file, f'{args.outdir}/warnings.txt', config)

    # select the best-matched sample sets
    logging.info('step2: select the best-matched sample sets')
    print(f'python3 {bindir}/deal_PoN.py -d {sample_depth_file} -o {args.outdir} -p {args.panel} -m {args.method} -c {args.config}\n')
    PoN_sets = PoN(sample_depth_file, PoN_list, args.outdir)
    PoN_sets.cal_delta() if args.method == 'delta' else PoN_sets.cal_corr()
    PoN_sets.select_ref_sets(config, ref_method=args.select_method)
    PoN_sets.construct_ref_matrix()

    # calculate ratio / z-score
    logging.info('step3: infer cnv')
    print(f'python3 {bindir}/infer_CNV.py -i {sample_depth_file} -r {refmatrix} -a {anno_file} -af {anno_file_cutoff} {gc_correct}\n')
    inferCNV = CNV(sample_depth_file, refmatrix, anno_file, anno_file_cutoff, config, args.gc_correct)
    inferCNV.do_CNV()
    inferCNV.filter_anno()
    print(f'grep -vP "chr13\t32973298|chr17\t41196282" {anno_file} > {anno_file_filter_UTR}') ## 删除UTR bed
    subprocess.call(f'grep -vP "chr13\t32973298|chr17\t41196282" {anno_file} > {anno_file_filter_UTR}', shell=True)

    # merge cnv
    logging.info('step4: join adjacent exon')
    print(f'python3 {bindir}/merge_CNV.py -a {anno_file_filter_UTR} -o {args.outdir}\n')
    do_mergeCNV = mergeCNV(anno_file_filter_UTR, args.outdir, config)
    do_mergeCNV.do_merge()

    # select brca1/2 result
    logging.info('step5: get final LGR')
    print(f'python3 {bindir}/get_Report.py -i {exon_level_cnv} -r {args.refGeneWithVer} -o {report_file}\n')
    transcript_dict = getExon_num(args.refGeneWithVer)
    getReport(exon_level_cnv, report_file, transcript_dict)

    # plot fig
    logging.info('step6: plot fig')
    print(f'{config.get("software", "Rscript")} {bindir}/plot_102_copy_num-v2.R {anno_file} {fig_cn_pdf} && {config.get("software", "convert")} -quality 100 -density 360 {fig_cn_pdf} {fig_cn_png}')
    subprocess.call(f'{config.get("software", "Rscript")} {bindir}/plot_102_copy_num-v2.R {anno_file} {fig_cn_pdf} && {config.get("software", "convert")} -quality 100 -density 360 {fig_cn_pdf} {fig_cn_png}', shell=True)
    
if __name__ == "__main__":
    main()

