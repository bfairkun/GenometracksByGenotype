#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
#
#
######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : SumBigwigs
# @created     : Monday Oct 17, 2022 13:01:12 CDT
#
# @description : 
######################################################################

import sys
import numpy
import pyBigWig



# I like to script and debug with an interactive interpreter. If using
# interactive interpreter to script quick args with sys.argv, parse_args with
# hardcoded args below


def SummariseBigwigsInWindow(bigwig_out_fn, bigwigs_in_fn_list, window, function="sum"):
    """
    Given output bigwig filepath, and list of input bigwig filepaths, summarise the input bigwigs with a function. For example, sum the values in all the input files. bigwig chrom headers should match.
    Note that this reads all bigwigs_in_fn_list file data into memory by iteration_windowsize_bp, so for memory and computational efficiency, consider changing the windowsize according to sparsity of data. This is best done for relatively sparse bigwigs (only data for a small number of loci)
    """
    chrom, coords = window.split(':')
    coords_float = [int(i.replace(',', '')) for i in coords.split('-')]
    print(coords_float)
    # check that input bigwig headers chrom list matches across all files
    bigwigs = [pyBigWig.open(bw) for bw in bigwigs_in_fn_list]
    chroms = bigwigs[0].chroms()
    for i, bw in enumerate(bigwigs):
        if bw.chroms() != chroms:
            raise Exception(f'bigwig chromosome header list should match. {bigwigs_in_fn_list[0]} chroms do not match {bigwigs_in_fn_list[i]}.')
    # Initialize output file
    bw_out = pyBigWig.open(bigwig_out_fn, "w")
    bw_out.addHeader(list(bw.chroms().items()))
    # Iterate through chromosomes and write out summarise bigwig intervals
    Array = []
    for bw in bigwigs:
        if bw.stats(chrom, *coords_float, type="max") == [None]:
            continue
        else:
            Array.append(bw.values(chrom, *coords_float,numpy=True))
    if Array:
        print("Writing", chrom, *coords_float)
        ArraySummarised = numpy.sum(numpy.array(Array), axis=0)
        bw_out.addEntries(chrom, coords_float[0], values=ArraySummarised, span=1, step=1)   
    bw_out.close()


if __name__ == '__main__':
    if hasattr(sys, 'ps1'):
        sys.argv = ["", "scratch/test.bw" , "chr2:197,623,164-198,410,427", "bigwigs/H3K4ME3/NA18489.1.bw", "bigwigs/H3K4ME3/NA18498.2.bw"]

    _, out_fn, Region = sys.argv[0:3]
    bigwigs_in = sys.argv[3:]
    SummariseBigwigsInWindow(out_fn, bigwigs_in, Region)

