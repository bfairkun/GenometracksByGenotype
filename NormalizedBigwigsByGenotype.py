"""
Outputs:
    <OutputPrefix>_output_0.bw
    <OutputPrefix>_output_1.bw
    <OutputPrefix>_output_2.bw
    <OutputPrefix>_output_tracks.ini
    <OutputPrefix>_PerInd_<Genotype>_<Sample>.bw

Note that bigwigs must not contain NaN values. Use 0 to denote no coverage. This can be accomplished with bedtools genomecov with -bga option then convert bedGraphToBigWig. See helper script in this repo.
--File or autodetect for samples to files
"""

import sys
import os
import argparse
import pandas as pd
import vcf
import pyBigWig
import numpy
import glob
import re
from jinja2 import Template
from collections import defaultdict

def ReadSamplesTsv(fn):
    """
    read tab delimited file of sample ID (must match vcf sample if attempting
    to average by genotype), bigwig filepath, group_label (optional), and 0/1
    for non-inverted/inverted (optional, default is non-inverted) which can be
    useful for plotting stranded coverage
    """
    return()


def str_to_int(MyStr):
    return(int(MyStr.replace(",", "")))

def GetNormalizationFactorBasedOnRegionList(bigwig_fn, regionlist):
    """
    regionlist is a list of lists. Each list entries is formatted like this:
    ["chr1", 1, 10].
    """
    try:
        bw = pyBigWig.open(bigwig_fn)
        TotalSignal = []
        for entry in regionlist:
            EntrySignal = sum(bw.values(RegionChr, RegionStart, RegionStop, numpy=True))
            TotalSignal.append(EntrySignal)
        return(sum(TotalSignal)/1E3)
    except:
        return(None)

def invert_dict(d): 
    inverse = dict() 
    for key in d: 
        val = d[key] 
        if val not in inverse: 
            inverse[val] = [key] 
        else: 
            inverse[val].append(key) 
    return inverse

def WriteLinksFileByGenotypes(SampleToGenotypeDict, bed_df, f_out_prefix=""):
    GenotypesToSamplesDict = invert_dict(SampleToGenotypeDict)
    bed_df.rename(columns={ bed_df.columns[0]: "#Chr", bed_df.columns[1]:"start", bed_df.columns[2]: "end", bed_df.columns[3]: "pid" }, inplace = True)
    for genotype in GenotypesToSamplesDict.keys():
        df_out = bed_df[bed_df.columns.intersection(GenotypesToSamplesDict[genotype] + ['#Chr','start', 'end', 'pid'])]
        df_out['psi'] = df_out.drop(['#Chr', 'start', 'end', 'pid'], axis=1).mean(axis=1).apply(lambda x:"{0:.2g}".format(x) )
        df_out.loc[:, ['#Chr',  'start', 'end', '#Chr', 'start', 'end','psi']].to_csv(f_out_prefix + "output_" + str(genotype) + ".links", sep='\t', index=False, header=False)

def GetNormalizationFactorWholeGenome(bigwig_fn):
    try:
        bw = pyBigWig.open(bigwig_fn)
        return(bw.header()['sumData']/1E9)
    except:
        return(None)

def parse_args(Args=None):
    p = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("VCF",
            help= "gzipped and tbi indexed vcf file with sample genotypes",
            metavar="<VCFFile>")
    p.add_argument("SnpPos",
            help= "Snp position. Position should be one-based, same as in vcf",
            metavar="<CHR>:<POS>")
    p.add_argument("Region",
            help= "Region to output",
            metavar="<CHR>:<START>-<STOP>")
    p.add_argument("BigwigList",
            help= """
            If using with the --BigwigListType KeyFile  option, a tab delimited text file with samples in the first column and a path to the input bigwig files in the second column. Alternatively, if using with the --BigwigListType GlobPattern option, use a  wildcard glob of bigwig files, and the samples will be inferred by searching the expanded filepaths for strings that match samples in the VCF. Note that using this wildcard method requires you to enclose the grob pattern in quotes, since this script needs to parse the grob pattern with python. Excluding quotes will immediately expand the grob pattern by the shell.

            Example1:
            MySampleToBigwigKey.tsv --BigwigListType KeyFile

            Example2:
            "./Bigwigs/*/Coverage.bw" --BigwigListType GlobPattern

            where the * expansion may be sample names, which are automatically matched (with some wiggle room). Using a keyfile is obviously more robust, as the glob pattern requires some guessing to match files to samples. For convenient use with 1000Genomes project samples, this sripts will match samples by searching filenames for 5 digit patterns that match sample naming scheme used in 1000Genomes project. For example: a file named ./MyBigwigs/GM19137.bw will be matched to the sample HG19137 because of the 19137 substring match""",
            metavar='{"<GlobPattern>",<KEYFILE>}')
    p.add_argument('--BigwigListType',
                    choices=["KeyFile", "GlobPattern"],
                    help='Required. Define how the bigwig list positional argument is specified.', required=True)
    p.add_argument("--Normalization",
            help="A bed file of regions to use for sample-depth normalization. The within-sample total coverage over the regions will be used to normalize each sample. If 'None' is chosen, will not perform any normalization. If 'WholeGenome' is used, will use whole genome.  (default: %(default)s)",
            default='WholeGenome', metavar="{<BEDFILE>,None,WholeGenome}")
    p.add_argument("--OutputPrefix",
            help="Prefix for all output files (default: %(default)s)",
            default="./")
    p.add_argument("--BedfileForSashimiLinks",
            help="QTLtools style bed or bed.gz file with header for samples (sample names must match vcf) and values to calculate average PSI per genotype. (default: %(default)s)",
            default=None)
    p.add_argument("--OutputNormalizedBigwigsPerSample",
            help="Output normalized bigwigs for each sample.",
            action='store_true')
    p.add_argument("--TracksTemplate", metavar="<FILE>", help="A jinja template for a tracks file for pyGenomeTracks customization. An example is included. Template variables allowed are 'OutputPrefix', 'HomoRefTitle', 'HetTitle', 'HomoAltTitle', and 'YMax'. If this argument is provided, the template file will be populated with the template variables to create a tracks file that can be used for pyGenomeTracks. If this argument is not provided, will output a very basic tracks file that can be used for pyGenomeTracks")
    p.add_argument("-v", "--verbosity", type=int, choices=[0,1,2], default=0,
            help="increase output verbosity (default: %(default)s)")
    return(p.parse_args(Args))


def main(args):
    """
    main function of script. args is a argparse.Namespace object created by the parse_args function defined above. example usages:
    # run main with parameters defined from sys.argv
    main(parse_args())
    # run main with parameters defined in custom list:
    main(parse_args(['test_data/1KG.Subset.vcf.gz', 'chr4:173807', 'chr4:117320-138022', '"test_data/bigwigs/*"', '--Normalization', 'None', '--BigwigListType', 'GlobPattern', '--OutputPrefix', 'test_results/', '--TracksTemplate', 'tracks_templates/tracks.ini.template3.txt', '--OutputNormalizedBigwigsPerSample']))
    """
    SnpChr, SnpCoord = args.SnpPos.split(":")
    vcf_reader = vcf.Reader(filename=args.VCF)
    SnpList = []
    for record in vcf_reader.fetch(SnpChr, str_to_int(SnpCoord)-1, str_to_int(SnpCoord)):
        SnpList.append(record)
    if len(SnpList) <1:
        print("No SNP found at SnpPos")
    elif len(SnpList)>1:
        print("More than one SNP found at SnpPos. Make sure vcf has one line per position")
    else:
        RefAllele = SnpList[0].REF
        AltAllele = SnpList[0].ALT[0]
        print("Found SNP")

    # Save sample names, grouped by genotype, to lists
    hets = [call.sample for call in SnpList[0].get_hets()]
    homo_refs = [call.sample for call in SnpList[0].get_hom_refs()]
    homo_alts = [call.sample for call in SnpList[0].get_hom_alts()]

    # Dict to lookup genotype for each sample. Use 0,1,2 for shorthand coding for homoRef, het, homoAlt
    SampleToGenotypeDict = dict(zip(
        homo_refs + hets +  homo_alts,
        [0 for i in homo_refs] + [1 for i in hets] + [2 for i in homo_alts]))
    # print(SampleToGenotypeDict)

    if args.BedfileForSashimiLinks:
        bed = pd.read_csv(args.BedfileForSashimiLinks, sep='\t')
        WriteLinksFileByGenotypes(SampleToGenotypeDict, bed, f_out_prefix=args.OutputPrefix)


    #If key file, make BigwigToSampleDict
    #Check that each file exists, warn if not

    # Get bigwig grob pattern, try to match to samples
    BigwigGlobmatchs = glob.glob(args.BigwigList)
    AllGenotypedNumeric = [re.search("\d{5}", sample).group(0) for sample in SampleToGenotypeDict.keys()]
    NumericToSamplenameDict = dict(zip(AllGenotypedNumeric, SampleToGenotypeDict.keys()))

    BigwigToSampleDict = defaultdict(dict)
    for fn in BigwigGlobmatchs:
        Match = re.findall("\d{5}", fn)
        if len(Match) == 1 and Match[0] in AllGenotypedNumeric:
            # BigwigToNumericDict[fn] = Match[0]
            # print(Match[0])
            SampleName = NumericToSamplenameDict[Match[0]]
            BigwigToSampleDict[fn]["SampleName"]= SampleName
            BigwigToSampleDict[fn]["Genotype"]=SampleToGenotypeDict[SampleName]


    #TODO: write out inferred sampleKeyFile

    #Check if bedfile given for normalization, and read regions if so.
    if args.Normalization not in ['None', 'WholeGenome']:
        try:
            NormalizationRegions = []
            with open(args.Normalization, 'r') as fh:
                for line in fh:
                    l = line.strip('\n').split('\t')
                    NormalizationRegions.append(l[0:3])
        except:
            print("Normalization bedfile non-existent or corrupted")

    #Iterate through bigwigs, add normalized values row by row to 2D array (list of 1D arrays) for each genotype.
    RegionChr, RegionCoords = args.Region.split(":")
    RegionStart, RegionStop = [str_to_int(i) for i in RegionCoords.split("-")]
    DictOfArraysForEachGenotype={0:[], 1:[], 2:[]}
    # print(BigwigToSampleDict)
    for i,bigwig in enumerate(BigwigToSampleDict.keys()):
        # Get normalization factor for each bigwig
        try:
            if args.Normalization == 'WholeGenome':
                NormFactor = GetNormalizationFactorWholeGenome(bigwig)
            elif args.Normalization == 'None':
                NormFactor = 1
            else:
                NormFactor = GetNormalizationFactorBasedOnRegionList(bigwig, NormalizationRegions)
        except:
            NormFactor = None
        if NormFactor:
            genotype = BigwigToSampleDict[bigwig]['Genotype']
            bw = pyBigWig.open(bigwig)
            NormalizedValues = bw.values(RegionChr, RegionStart, RegionStop, numpy=True)/NormFactor
            DictOfArraysForEachGenotype[genotype].append(NormalizedValues)
            bwOut_fn = "output_PerInd_{}_{}.bw".format(genotype, BigwigToSampleDict[bigwig]["SampleName"])
            BigwigToSampleDict[bigwig]["bwOut_fn"] = bwOut_fn
            if args.OutputNormalizedBigwigsPerSample:
                bwOut = pyBigWig.open(args.OutputPrefix + bwOut_fn, "w")
                bwOut.addHeader(list(bw.chroms().items()))
                bwOut.addEntries(RegionChr, RegionStart, values=NormalizedValues, span=1, step=1)
                bwOut.close()

    # average the normalized signal across samples by genotype and write out bigwigs and tracks file
    MaxValuesOut = []
    MaxValuesPerInd = []
    for genotype, Array in DictOfArraysForEachGenotype.items():
        if len(Array) > 0:
            ArrayOut = numpy.mean(numpy.array(Array), axis=0)
            bwOut = pyBigWig.open(args.OutputPrefix + "output_" + str(genotype) + "_.bw", "w")
            bwOut.addHeader(list(bw.chroms().items()))
            bwOut.addEntries(RegionChr, RegionStart, values=ArrayOut, span=1, step=1)
            bwOut.close()
            #Save the max signal in window for to set ylimits in tracks file
            MaxValuesOut.append(max(ArrayOut))
            MaxValuesPerInd.append(numpy.max(Array))
        else:
            print("Genotype not found")
    GenotypeCounts = {k:len(v) for (k, v) in DictOfArraysForEachGenotype.items()}
    YMax = max(MaxValuesOut)
    YMax_PerInd = max(MaxValuesPerInd)


    #Get ind output bigwigs separated by genotype for access as jinja variable
    HomoRefBwList = []
    HetBwList = []
    HomoAltBwList = []
    for k,v in BigwigToSampleDict.items():
        try:
            bigwig = v['bwOut_fn']
        except KeyError:
            continue
        if v['Genotype'] == 0:
            HomoRefBwList.append(bigwig)
        elif v['Genotype'] == 1:
            HetBwList.append(bigwig)
        elif v['Genotype'] == 2:
            HomoAltBwList.append(bigwig)

    if args.TracksTemplate:
        with open(args.TracksTemplate) as file_:
            template_file_contents = file_.read()
            template = Template(template_file_contents)
    else:
        if args.BedfileForSashimiLinks:
            template = Template(
    """
    [output_0_]
    bw_file = {{OutputPrefix}}output_0_.bw
    link_file = {{OutputPrefix}}output_0.links
    title = {{HomoRefTitle}}
    height = 2
    bw_color = red
    min_value = 0
    max_value = {{YMax}}
    number_of_bins = 2000
    nans_to_zeros = true
    show_data_range = true
    y_axis_values = original
    link_color = black
    line_style = solid
    fontsize = 8
    # The link in Sashimi plot is a Bezier curve.
    # The height of the curve is calculated from the length of the intron.
    # When the y-axis in bigwig track is different, the height of curve needs to be scaled.
    scale_link_height = 1
    # The line width for links is proportion to the numbers at the last column in links file (PSI).
    # But the absolute width is calculated from the supplied numbers, which can look too thin or too wide sometimes.
    # Use scale_line_width to scale the absolute line widths.
    # You may need to try several values to get a satisfying result.
    scale_line_width = 0.05
    file_type = sashimiBigWig
    show_number = true

    [output_1_]
    bw_file = {{OutputPrefix}}output_1_.bw
    link_file = {{OutputPrefix}}output_1.links
    title = {{HetTitle}}
    height = 2
    bw_color = purple
    min_value = 0
    max_value = {{YMax}}
    number_of_bins = 700
    nans_to_zeros = true
    show_data_range = true
    y_axis_values = original
    link_color = black
    line_style = solid
    fontsize = 8
    scale_link_height = 1
    scale_line_width = 0.05
    file_type = sashimiBigWig
    show_number = true

    [output_2_]
    bw_file = {{OutputPrefix}}output_2_.bw
    link_file = {{OutputPrefix}}output_2.links
    title = {{HomoAltTitle}}
    height = 2
    bw_color = blue
    min_value = 0
    max_value = {{YMax}}
    number_of_bins = 700
    nans_to_zeros = true
    show_data_range = true
    y_axis_values = original
    link_color = black
    line_style = solid
    fontsize = 8
    scale_link_height = 1
    scale_line_width = 0.05
    file_type = sashimiBigWig
    show_number = true

    [vlines]
    type = vlines
    file = {{OutputPrefix}}output_SNP.bed
    """
    )
        else:
            template = Template(
    """
    [output_0_avg]
    file = {{OutputPrefix}}output_0_.bw
    height = 2
    color = red
    type = line
    min_value = 0
    max_value = {{YMax}}
    number_of_bins = 2000
    nans_to_zeros = true
    show_data_range = true
    y_axis_values = original
    file_type = bigwig

    [output_1_avg]
    file = {{OutputPrefix}}output_1_.bw
    height = 2
    color = purple
    type = line
    min_value = 0
    max_value = {{YMax}}
    number_of_bins = 2000
    nans_to_zeros = true
    show_data_range = true
    y_axis_values = original
    file_type = bigwig
    overlay_previous = share-y


    [output_2_avg]
    file = {{OutputPrefix}}output_2_.bw
    height = 2
    color = blue
    type = line
    min_value = 0
    max_value = {{YMax}}
    number_of_bins = 2000
    nans_to_zeros = true
    show_data_range = true
    y_axis_values = original
    file_type = bigwig
    overlay_previous = share-y


    [vlines]
    type = vlines
    file = {{OutputPrefix}}output_SNP.bed
    labels = true
    """
    )

    with open(args.OutputPrefix + "output_tracks.ini", 'w') as file_out:
        print("writing ini file")
        _ = file_out.write(
            template.render(
                        #Jinja variables
                        OutputPrefix=args.OutputPrefix,
                        YMax = YMax, #Max after averaging by genotype
                        YMax_PerInd = YMax_PerInd, #Max across all individuals
                        RefAllele = RefAllele,
                        AltAllele = AltAllele,
                        GenotypeCounts = GenotypeCounts, #Dictionary, with keys of 0,1,2
                        HomoRefBwList = HomoRefBwList, #List of bigwigs for each ind with homoRef genotype
                        HetBwList = HetBwList,
                        HomoAltBwList = HomoAltBwList,
                        HomoRefTitle="0: {Ref}/{Ref} (n={count})".format(Ref=RefAllele, count=GenotypeCounts[0]),
                        HetTitle="1: {Ref}/{Alt} (n={count})".format(Ref=RefAllele, Alt=AltAllele, count=GenotypeCounts[1]),
                        HomoAltTitle="2: {Alt}/{Alt} (n={count})".format(Alt=AltAllele, count=GenotypeCounts[2])))

    with open(args.OutputPrefix + "output_SNP.bed", "w") as file_snp_bed:
        _ = file_snp_bed.write("{}\t{}\t{}\t{}".format(SnpChr, str_to_int(SnpCoord)-1, SnpCoord, SnpList[0].ID))

if __name__ == '__main__':
    # I like to script and debug with an interactive interpreter. If using
    # interactive interpreter, can step thru the script with args defined below
    try:
        if(sys.ps1):
            Args = 'test_data/1KG.Subset.vcf.gz chr4:173807 chr4:117320-138022 "test_data/bigwigs/*" --Normalization None --BigwigListType GlobPattern --OutputPrefix test_results/ --TracksTemplate tracks_templates/tracks.ini.template3.txt --OutputNormalizedBigwigsPerSample'.split(' ')
            args = parse_args(Args=Args)
    except:
        args = parse_args()
    main(args)

