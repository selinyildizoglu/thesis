#!/usr/bin/env python
# encoding: utf-8

# Reads BAM files and extracts alleles at specific positions

import pandas as pd
import pysam

from collections import defaultdict as ddict
import sys, os, glob


def print_log(msg, verbose=True):
    if verbose:
        print(msg)

def filter_snps(snps, verbose=True):
    pos = []
    for x in snps.fetch():
        alt = x.alts[0]
        if not (len(x.ref) < 2) or not (len(alt) < 2):
            print_log(
                f'> Position {x.pos} discarded - structural variant',
                verbose)
            continue
        dp4 = x.info['DP4']
        ref_prop = (dp4[0] + dp4[1])/sum(dp4)
        pos_dict = {'pos' : x.pos, 'ref' : x.ref,
                'alt' : alt,'ref_prop_exp' : ref_prop}
        pos.append(pos_dict)
    return pos

def extract_bases(alignments, positions, ref='BK000964.3_looped_3008'):
    bases_reads = []
    for pos in positions:
        for pcol in alignments.pileup(ref, pos-1, pos):
            if pcol.pos != (pos-1):
                continue
            for pread in pcol.pileups:
                if not pread.is_del and not pread.is_refskip:
                    bases_reads.append({
                        'read_id' : pread.alignment.query_name,
                        'snp_pos' : pos,
                        'snp_allele' : pread.alignment.query_sequence[
                            pread.query_position]})
    return bases_reads

def select_reads(snp_calls, meth_calls):
    return snp_calls[snp_calls.read_id.isin(
        meth_calls.read_id.values)].drop_duplicates()

def _calculate_coverage(counts):
    return counts.sum(axis=1)

def _find_covered_positions(counts, min_cov=10):
    cov = _calculate_coverage(counts)
    return counts.loc[cov >= min_cov].index.values

def _get_alleles_pos(base_counts, maf_th=0.1):
    return [(allele, base_counts[allele])
            for allele in base_counts.index
            if base_counts[allele] >= maf_th*sum(base_counts)]

def _get_alleles(counts, maf_th=0.1):
    alleles = { pos : _get_alleles_pos(base_counts, maf_th=maf_th)
            for pos, base_counts in counts.iterrows()}
    return alleles


def select_ref(x):
    return x[x['ref']]

def _parse_args():
    import argparse

    parser = argparse.ArgumentParser(
        description='SNP Alleles.')
    parser.add_argument('sample', metavar='S', type=str,
                               help='root name of sample to process')
    parser.add_argument('-r', '--region', type=str, default="BK000964.3_looped_3008",
                               help='contig name')
    parser.add_argument('-o', '--output', type=str,
                               help='path to output folder')
    parser.add_argument('-a', '--align', type=str,
                               help='path to alignments folder')
    parser.add_argument('-s', '--snps', type=str,
                               help='path to SNPs folder')
    parser.add_argument('-c', '--coverage', type=int,
                               help='minimum coverage')
    parser.add_argument('-t', '--threshold', type=float,
                               help='minimum allele frequency threshold')
    parser.add_argument('-v', '--verbose', action='store_true',
                               help='increase output verbosity')

    return parser.parse_args()

def main():


# 00.- Configuration

    args=_parse_args()
    root_fn = args.sample
    region=args.region

    print(f"Starting execution for sample {root_fn}")

## Paths

    ofd = f"{args.output}/{root_fn}/"

    if not os.path.exists(ofd):
            os.makedirs(ofd)

    als_fd = f"{args.align}/{root_fn}/"
    als_fn = glob.glob(f"{als_fd}*.bam")[0]

    if os.path.isdir(args.snps):
        snps_fd = f"{args.snps}/{root_fn}/"
        # TO DO: Remove this
        snps_fn = glob.glob(f"{snps_fd}*.strain.vcf")[0]
    else:
        snps_fn = args.snps

## Options

    min_cov = args.coverage
    maf_th = args.threshold
    verbose = args.verbose

# 01.- Read SNPs

    snps = pd.DataFrame(filter_snps(pysam.VariantFile(snps_fn, 'r')))

# 02.- Read alignment files

    als = pysam.AlignmentFile(als_fn, 'rb')

# 03.- Extract bases from reads

    calls = pd.DataFrame(
            extract_bases(als, snps['pos'], ref=region)).sort_values(
            by=['read_id', 'snp_pos'])

    calls.to_csv(f"{ofd}{root_fn}_alleles.csv",
            index=False)

if __name__ == "__main__":
    main()
