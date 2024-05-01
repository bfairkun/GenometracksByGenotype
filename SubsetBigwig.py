#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
#
#
######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : SubsetBigwig
# @created     : Tuesday Apr 30, 2024 13:18:48 CDT
#
# @description : 
######################################################################

import pyBigWig
import argparse

def subset_bigwig(input_bigwig, output_bigwig, region):
    # Parse the region string
    chrom, region_range = region.split(':')
    start, end = region_range.split('-')
    start, end = int(start.replace(',', '')), int(end.replace(',', ''))

    # Open the input BigWig file
    bw_in = pyBigWig.open(input_bigwig)

    # Open a new BigWig file for writing
    bw_out = pyBigWig.open(output_bigwig, 'w')

    # Get the chromosome sizes from the input BigWig file
    chrom_sizes = [(c, bw_in.chroms()[c]) for c in bw_in.chroms()]

    # Check if the specified chromosome exists
    if chrom not in bw_in.chroms():
        raise ValueError("Chromosome {} not found in the input BigWig file.".format(chrom))

    # Subset the data
    subset_values = bw_in.values(chrom, start, end)

    # Prepare lists of chromosomes, starts, ends, and values
    chrom_list = [chrom] * len(subset_values)
    start_list = list(range(start, end))
    end_list = list(range(start + 1, end + 1))

    # Write subset to the output BigWig file
    bw_out.addHeader(chrom_sizes)
    bw_out.addEntries(chrom_list, start_list, end_list, values=subset_values)


    # Close the BigWig files
    bw_in.close()
    bw_out.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Subset a BigWig file to a specified genomic region")
    parser.add_argument("input_bigwig", help="Input BigWig file path")
    parser.add_argument("output_bigwig", help="Output BigWig file path")
    parser.add_argument("region", help="Genomic region to subset (format: chr:start-end)")
    args = parser.parse_args()

    subset_bigwig(args.input_bigwig, args.output_bigwig, args.region)
