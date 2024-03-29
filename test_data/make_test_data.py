#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
#
#
######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : make_test_data
# @created     : Wednesday Aug 18, 2021 12:48:40 CDT
#
# @description :
######################################################################

import sys
# I like to script and debug with an interactive interpreter. If using
# interactive interpreter to script quick args with sys.argv, parse_args with
# hardcoded args below
if hasattr(sys, 'ps1'):
    sys.argv = ["", "Arg1" ,"Arg2"]

import pyBigWig
import numpy
import os

bw_list = os.listdir("/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/bigwigs/Expression.Splicing/")

DirOut = "test_bigwigs/"
os.makedirs(DirOut, exist_ok=True)

region = "chr4:39,031,787-39,110,467"
chrom,coords = region.split(":")
start, stop = [int(i) for i in coords.replace(",", "").split("-")]

basename = bw_list[0]

for i, basename in enumerate(bw_list):
    # if i == 0:
    bw_out = pyBigWig.open(DirOut + "/" +  basename, "w")
    bw = pyBigWig.open("/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/bigwigs/Expression.Splicing/" + basename)
    bw_out.addHeader(list(bw.chroms().items()))
    for x in bw.intervals(chrom, start, stop):
        bw_out.addEntries([chrom],[x[0]],ends=[x[1]],values=[x[2]])
    bw.close()
    bw_out.close()


