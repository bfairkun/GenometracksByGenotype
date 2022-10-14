
import sys
import os
import tempfile
import glob
import logging
import argparse

import pandas as pd
import numpy
import pyBigWig

import AggregateBigwigsForPlotting


def SplitNIntoWindows(N, windowsize):
    """
    If N=20 and windowsize=7,
    return [(0,7), (7,14), (14,20)]
    """
    range_list = list(range(0, N, windowsize))
    if range_list[-1] != N:
        range_list.append(N)
    OutputList = []
    for i, n in enumerate(range_list[0:-1]):
        OutputList.append((n, range_list[i+1]))
    return(OutputList)

def SummariseBigwigs(bigwig_out_fn, bigwigs_in_fn_list, function="sum", iteration_windowsize_bp=1000000):
    """
    Given output bigwig filepath, and list of input bigwig filepaths, summarise the input bigwigs with a function. For example, sum the values in all the input files. bigwig chrom headers should match.
    Note that this reads all bigwigs_in_fn_list file data into memory by iteration_windowsize_bp, so for memory and computational efficiency, consider changing the windowsize according to sparsity of data. This is best done for relatively sparse bigwigs (only data for a small number of loci)
    """
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
    for ChromName, ChromLen in chroms.items():
        Windows = SplitNIntoWindows(ChromLen, iteration_windowsize_bp)
        # Iterate by windowsize
        for window in Windows:
            Array = []
            for bw in bigwigs:
                if bw.stats(ChromName, *window, type="max") == [None]:
                    continue
                else:
                    Array.append(bw.values(ChromName, *window, numpy=True))
            if Array:
                print("Writing", ChromName, window)
                ArraySummarised = numpy.sum(numpy.array(Array), axis=0)
                bw_out.addEntries(ChromName, window[0], values=ArraySummarised, span=1, step=1)   
    bw_out.close()

def parse_args(Args=None):
    p = argparse.ArgumentParser(
        conflict_handler="resolve",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
    "--QTLsBed",
    help="""bed file containing regions to plot. Must contain 6 columns. chrom, start, stop, SnpPos (ie chr7:1234), signed float (where positive sign indicates the REF is low allele), strand""",
    metavar='<BEDFILE>',
    )
    p.add_argument(
        "--BigwigList",
        help="""a tab delimited text file with samples in the first column and a path to the input bigwig files in the second column.""",
        metavar='{"<GlobPattern>",<KEYFILE>}',
    )
    p.add_argument(
        "--VCF",
        help="gzipped and tbi indexed vcf file with sample genotypes",
        metavar="<VCFFile>",
        default=None,
    )
    p.add_argument(
        "--Normalization",
        help="A bed file of regions to use for sample-depth normalization. The within-sample total coverage over the regions will be used to normalize each sample. If 'None' is chosen, will not perform any normalization. If 'WholeGenome' is used, will use whole genome. I haven't ever tested the <BEDFILE> option. (default: %(default)s)",
        default="None",
        metavar="{<BEDFILE>,None,WholeGenome}",
    )
    p.add_argument(
        "--FlankingRegionLengthBp",
        help="Number of base pairs to add to each side of region to include in output bigwigs. (default: %(default)s)",
        default="0",
        type=int,
        metavar="<Integer>",
    )
    p.add_argument(
        "--OutputPrefix",
        help="Prefix for all output files (default: %(default)s)",
        default="./",
    )
    p.add_argument(
        "--GroupSettingsFile",
        help=""" Requires use of --BigwigListType KeyFile. This option specifies tab delimited file with one row per Group_label to add extra plotting settings for each group. Columns are as follows: 1: Group_label (must match a Group_label in the KeyFile) 2: Group_color (optional). Hex or rgb colors to plot potentially plot each group as defined in the output ini 3: BedgzFile (optional). Any bedgz file provided in this column will take precedent over bedgz files provided to the --BedfileForSashimiLinks option. If a group settings file is provided only groups described in that file will be plotted. Additionally, tracks will be plotted in the order they are listed in this file""",
        default=None,
    )
    p.add_argument(
        "--Workdir",
        metavar="<path>",
        help="An optional path to set as workdir before executing main script function. Could be useful, for example, if bigwig file list uses relative file paths from a different directory",
        default="./",
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

def main(args):
    """
    return dataframe of samples, bigwig files, and do main work, of this script, resulting in output files including summarised bigwigs for each group and GenotypeByEffect (High, Med, Low).
    """
    logging.debug(args)
    DF_QTLs = pd.read_csv(args.QTLsBed, sep='\t', header=0, names=["chrom", "start", "stop", "SNPPos", "signedEffect", "strand"])
    DF_BigwigsForEachQTL = pd.DataFrame()
    temp_dir = tempfile.TemporaryDirectory()
    logging.debug(temp_dir.name)
    # execute AggregateBigwigsForPlotting to write bigwigs
    # And append to lists of bigwig files that will need to be merged to create a PosPos.bw, PosNeg.bw, and NegNeg.bw 
    for index, row in DF_QTLs.iterrows():
        print(row)
        MyArgs = f"--Workdir {args.Workdir} --Region {row['chrom']}:{row['start']-args.FlankingRegionLengthBp}-{row['stop']+args.FlankingRegionLengthBp} --SnpPos {row['SNPPos']} --VCF {args.VCF} --BigwigListType KeyFile --GroupSettingsFile {args.GroupSettingsFile} --BigwigList {args.BigwigList} --OutputPrefix {temp_dir.name}/{row['chrom']}-{row['start']}-{row['stop']}.  --NoSashimi -vv --TracksTemplate /project2/yangili1/bjf79/GenometracksByGenotype/tracks_templates/GeneralPurposeColoredByGenotype.ini"
        print(MyArgs)
        parsed_args = AggregateBigwigsForPlotting.parse_args(filter(None, MyArgs.split(' ')))
        DF_temp = AggregateBigwigsForPlotting.main(**vars(parsed_args))
        DF_temp["SNPPos"]=row['SNPPos']
        DF_temp["Region"]= f"{row['chrom']}_{row['start']}_{row['stop']}"
        DF_temp["signedEffect"]=row['signedEffect']
        if row['signedEffect'] > 0:
            DF_temp['GenotypeDefinedByEffect'] = DF_temp['genotype'].map({'0':"High", '1':"Mid", '2':"Low", '3':"Ungenotyped"})
        elif row['signedEffect'] < 0:
            DF_temp['GenotypeDefinedByEffect'] = DF_temp['genotype'].map({'0':"Low", '1':"Mid", '2':"High", '3':"Ungenotyped"})
        DF_BigwigsForEachQTL = pd.concat([DF_BigwigsForEachQTL, DF_temp])
    # Then merge the bigwigs by GenotypeDefinedByEffect (which is determined by the sign of the signedEffect) and write out new bigwig.
    DF_BigwigsForEachQTL_grouped = DF_BigwigsForEachQTL[["Group_label", "Strand","GenotypeDefinedByEffect", "GroupsSummarisedBigwigOut"]].groupby(["Group_label", "Strand","GenotypeDefinedByEffect"], as_index=False)[["GroupsSummarisedBigwigOut"]].agg(lambda x: list(x))
    for index, row in DF_BigwigsForEachQTL_grouped.iterrows():
        f_out = f'{args.OutputPrefix}Summarised-{row["Group_label"]}-{row["Strand"]}-{row["GenotypeDefinedByEffect"]}.bw'
        print(f"Writing {f_out}...")
        SummariseBigwigs(f_out, row["GroupsSummarisedBigwigOut"])
    # rm temp_dir when done:
    temp_dir.cleanup()
    return DF_BigwigsForEachQTL

if __name__ == "__main__":
    # I like to script and debug with an interactive interpreter in the same
    # working dir as this script.. If using interactive interpreter, can step
    # thru the script with args defined below
    try:
        if sys.ps1:
            Args = "--QTLsBed /project2/yangili1/bjf79/ChromatinSplicingQTLs/code/scratch/TestMetaplot.bed --FlankingRegionLengthBp 10000 --Workdir /project2/yangili1/bjf79/ChromatinSplicingQTLs/code/ --VCF /project2/yangili1/bjf79/ChromatinSplicingQTLs/code/Genotypes/1KG_GRCh38/Autosomes.vcf.gz --OutputPrefix /project2/yangili1/bjf79/ChromatinSplicingQTLs/code/scratch/TestMetaplotDir/ --BigwigList /project2/yangili1/bjf79/ChromatinSplicingQTLs/code/scratch/bwList.tsv --GroupSettingsFile /project2/yangili1/bjf79/ChromatinSplicingQTLs/code/scratch/TestMetaplotDir/Groups.tsv".split(" ")
            args = parse_args(Args=Args)
    except:
        args = parse_args()
    AggregateBigwigsForPlotting.setup_logging(args.verbosity)
    DF = main(args)
