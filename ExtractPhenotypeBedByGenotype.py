"""
Input: bedtools style phenotype table, and vcf, and snp, and (optional region, and feature name).
Outputs: a table, one row per sample, each column as different features, including the genotype at the specified SNP. Features output are the set intersection of all features within the same chromosome as the required SnpPos argument, that are within the optional 'CisWindowAroundSnpToGetFeatures' argument, and whose 4th column in the Bed file matches the optional 'FeatureName' argument.
"""

import sys
import os
import tempfile
import glob
import logging
import argparse

import pandas as pd
import numpy

import pysam

import AggregateBigwigsForPlotting


def parse_args(Args=None):
    p = argparse.ArgumentParser(
        conflict_handler="resolve",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--Bed",
        help="QTLTools formatted bgzipped and indexed bedfile",
        metavar="<BedgzFile>",
        required=True,
    )
    p.add_argument(
        "--CisWindowAroundSnpToGetFeatures",
        help="Cis window, in bp, around SnpPos to lookup features. If none is provided, will use the intersection of all features on the same chromosome as SnpPos, and the feature provided in FeatureName argument. Providing a window may slightly decrease lookup time of features",
        metavar="<int>",
        required=False,
    )
    p.add_argument(
        "--VCF",
        help="gzipped and tbi indexed vcf file with sample genotypes",
        metavar="<VCFFile>",
        default=None,
    )
    p.add_argument(
        "--SnpPos",
        help="Snp position. Position should be one-based, same as in vcf",
        metavar="<CHR>:<POS>",
        default=None,
    )
    p.add_argument(
        "--SnpName",
        help="Snp name. Name should have the following format: chr:position:Ref:Alt",
        metavar="<CHR>:<POS>:<REF>:<ALT>",
        default=None,
        required=False,
    )
    p.add_argument(
        "--FeatureName",
        help="Feature name. Should match col4 of bed file for specific feature to select. Otherwise, all features in overlapping region will be used",
        metavar="<FeatureName>",
        default=None,
        required=False,
    )
    p.add_argument(
        "--Out",
        help="Output file for tsv. Stdout if non provided",
        metavar="<FileOut>",
        default=None,
        required=False,
    )
    p.add_argument(
        "--Long",
        action="store_true",
        default=False,
        help="Output long format (useful if outputting multiple phenotypes)",
    )
    p.add_argument(
        "--ReportBedAsColumn",
        action="store_true",
        default=False,
        help="Add a column for the name of the Bed input file. Useful for concatenating long format output of repeated script calls with different Bed inputs files",
    )
    p.add_argument(
        "--NoHeader",
        action="store_true",
        default=False,
        help="Do not output header in output",
    )
    p.add_argument(
        "-v",
        "--verbose",
        action="count",
        dest="verbosity",
        default=0,
        help="verbose output (repeat for increased verbosity)",
    )
    p.add_argument(
        "-q",
        "--quiet",
        action="store_const",
        const=-1,
        default=0,
        dest="verbosity",
        help="quiet output (show errors only)",
    )
    return p.parse_args(Args)

def StringJoin(x, sep=";"):
    return sep.join([str(i) for i in x])

def main(args):
    """
    return dataframe of samples, bigwig files, and do main work, of this script, resulting in output files including summarised bigwigs for each group and GenotypeByEffect (High, Med, Low).
    """
    logging.debug(args)
    SnpChrom, SnpPos = args.SnpPos.split(':')
    if args.CisWindowAroundSnpToGetFeatures:
        RegionStart = int(SnpPos) - int(args.CisWindowAroundSnpToGetFeatures)
        RegionStop = int(SnpPos) + int(args.CisWindowAroundSnpToGetFeatures)
    else:
        RegionStart, RegionStop = None, None
    DF = AggregateBigwigsForPlotting.PysamTabixToPandas_df(args.Bed, SnpChrom, RegionStart, RegionStop)
    if args.FeatureName:
        DF = DF.loc[DF.iloc[:,3]==args.FeatureName]
    DF["Feature"] = DF.iloc[:,0:6].apply(StringJoin, axis=1 )
    DF = DF.drop(DF.columns[0:6], axis=1)
    DF = DF.set_index("Feature").transpose()
    DF.insert(loc = 0,column = 'SampleID',value = DF.index)
    DF = AggregateBigwigsForPlotting.AppendGenotypeColumnToDf(DF, vcf_fn=args.VCF, SnpPos=args.SnpPos, SnpName=args.SnpName)
    if args.Long:
        DF = pd.melt(DF, id_vars=['SampleID','genotype', 'Ref', 'Alt', 'ID'])
    if args.ReportBedAsColumn:
        DF["Bed"] = args.Bed
    if not args.Out:
        args.Out = sys.stdout
    DF.to_csv(args.Out, sep='\t', index=False, header=not args.NoHeader)
    return DF
# DF = main(args)

if __name__ == "__main__":
    # I like to script and debug with an interactive interpreter in the same
    # working dir as this script.. If using interactive interpreter, can step
    # thru the script with args defined below
    try:
        if sys.ps1:
            Args = "--Bed /project2/yangili1/bjf79/ChromatinSplicingQTLs/code/QTLs/QTLTools/Expression.Splicing/OnlyFirstRepsUnstandardizedForColoc.sorted.qqnorm.bed.gz --VCF /project2/yangili1/bjf79/ChromatinSplicingQTLs/code/Genotypes/1KG_GRCh38/Autosomes.vcf.gz -vvv --SnpName 14:105172594:T:C --SnpPos chr14:105172594".split(" ")
            args = parse_args(Args=Args)
    except:
        args = parse_args()
    AggregateBigwigsForPlotting.setup_logging(args.verbosity)
    DF = main(args)
