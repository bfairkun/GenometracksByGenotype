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
import pysam
import pyBigWig
import numpy
import glob
import re
import colorbrewer
import logging
from jinja2 import Template
from collections import defaultdict
from io import StringIO
import pybedtools
import gzip
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.backends.backend_pdf import PdfPages


pd.options.mode.chained_assignment = None

# Add a colorbrewer set with red, purple, blue, and gray. Nice colors for plotting genotypes
colorbrewer.RdPrBlu = {
    3: [(240, 59, 32), (122, 1, 119), (8, 81, 156)],
    4: [(240, 59, 32), (122, 1, 119), (8, 81, 156), (189, 189, 189)],
}

colorbrewer_palletes = [i for i in colorbrewer.__dict__.keys() if type(getattr(colorbrewer, i)) is dict and i != '__builtins__']

def str_to_int(MyStr):
    return int(MyStr.replace(",", ""))


def rgb_to_hex(rgb):
    return "#%02x%02x%02x" % rgb

def hex_to_rgb(hexcolor):
    h = hexcolor.lstrip('#')
    return tuple(int(h[i:i+2], 16) for i in (0, 2, 4))


def setup_logging(verbosity):
    """
    with base_loglevel=30...
    verbosity =-1 => ERROR
    verbosity = 0 => WARNING
    verbosity = 1 => INFO
    verbosity = 2 => DEBUG
    """
    base_loglevel = 30
    verbosity = min(verbosity, 2)
    loglevel = base_loglevel - (verbosity * 10)
    logging.basicConfig(
        level=loglevel,
        format="[%(asctime)s] {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s",
    )


def test_logging():
    logging.debug("debug message")
    logging.info("info message")
    logging.warning("warning message")
    logging.error("error message")
    logging.critical("critical message")


def ReadSamplesTsv(fn):
    """
    read 2+ column tab delimited file where the first 5 columns defined as:
    1: sample ID (must match vcf sample if attempting
    to average by genotype)
    2: bigwig filepath
    3: group_label (optional)
    4: ""|"+"|"-" for strand. (Do not actually enter quotes.. optional, default is "" for unstranded. this is useful for plotting stranded coverage.
    Columns past these 4 will be ignored
    """
    df = pd.read_csv(
        fn,
        header=0,
        names=["SampleID", "BigwigFilepath", "Group_label", "Strand"],
        sep="\t",
        dtype={"Strand": str},
    )
    df.fillna({"Group_label": "", "Strand": ""}, inplace=True)
    return df


def AddExtraColumnsPerGroup(
    df, GroupToBedTsv=None, BedfileForAll=None, MergeStyle="inner"
):
    """
    GroupToBed tsv file needed, otherwise, utilize make one on the fly with
    blank/default values
    """
    if GroupToBedTsv:
        GroupToBed_df = pd.read_csv(
            GroupToBedTsv,
            header=0,
            names=["Group_label", "Group_color", "BedgzFilepath", "Supergroup"],
            sep="\t",
        )
    else:
        GroupToBed_df = pd.DataFrame({'Group_label': df['Group_label'].unique(), 'Group_color':'', 'BedgzFilepath':'', 'Supergroup':''})
    # Set Supergroup to Group_label where Supergroup is empty or NaN
    mask = GroupToBed_df['Supergroup'].isna() | (GroupToBed_df['Supergroup'] == '')
    GroupToBed_df.loc[mask, 'Supergroup'] = GroupToBed_df.loc[mask, 'Group_label']
    GroupToBed_df['PlotOrder'] = numpy.arange(len(GroupToBed_df))
    if BedfileForAll:
        GroupToBed_df.fillna(BedfileForAll, inplace=True)
    logging.debug(GroupToBed_df)
    df = df.merge(GroupToBed_df, how=MergeStyle, on="Group_label")
    return df


def ConstructSamplesTsvFromGlobPattern(GlobPatternString, vcf_fn, fn_out=None):
    """
    Return df same as ReadSamplesTsv, with sample ID and bigwig filepath
    inferred from glob pattern. The group_label will default to "" for all
    samples, and Strand will default to "" (unstranded). If fn_out is
    specified, write out the df
    """
    BigwigFile_VcfSample_Matchlist = []
    BigwigGlobmatchs = glob.glob(GlobPatternString)
    vcf_reader = pysam.VariantFile(vcf_fn)
    vcf_samples = list(vcf_reader.header.samples)
    for bigwig_fn in BigwigGlobmatchs:
        # Get query to search amongst vcf samples
        Match = re.findall(r"\d{5}", bigwig_fn)
        if len(Match) == 1:
            query = Match[0]
        else:
            continue
        vcf_sample_matches = [i for i in vcf_samples if query in i]
        if len(vcf_sample_matches) == 1:
            vcf_sample_match = vcf_sample_matches[0]
        else:
            continue
        BigwigFile_VcfSample_Matchlist.append([vcf_sample_match, bigwig_fn])
    df = pd.DataFrame(
        BigwigFile_VcfSample_Matchlist, columns=["SampleID", "BigwigFilepath"]
    )
    df["Group_label"] = ""
    df["Strand"] = ""
    if fn_out:
        pass
    return df


def label_genotype(row, homo_alts, hets, homo_refs):
    if row["SampleID"] in homo_alts:
        return "2"
    elif row["SampleID"] in hets:
        return "1"
    elif row["SampleID"] in homo_refs:
        return "0"
    else:
        return "3"

def get_vcf_record(vcf_reader, SnpPos, SnpName=None):
    SnpChr, SnpCoord = SnpPos.split(":")
    SnpList = []
    len_records = 0
    # pysam uses 0-based coordinates
    fetch_start = str_to_int(SnpCoord) - 1
    fetch_end = str_to_int(SnpCoord)
    for record in vcf_reader.fetch(SnpChr, fetch_start, fetch_end):
        len_records += 1
    for record in vcf_reader.fetch(SnpChr, fetch_start, fetch_end):
        if SnpName and (len_records > 1):
            SnpChr2, SnpCoord2, SnpRef, SnpAlt = SnpName.split(':')
            if (record.ref == SnpRef) and (str(record.alts[0]) == SnpAlt):
                SnpList.append(record)
            else:
                continue
        else:
            SnpList.append(record)
    return SnpList

def AppendGenotypeColumnToDf(df, vcf_fn=None, SnpPos=None, SnpName=None):
    """
    Append a column of 0, 1, or 2 for homozygous ref, het, or homozygous alt
    genotype for for each entry in df, based on the SampleID column, the
    vcf_fn, and the SNP location. 3 is returned in other cases (no genotype
    found)
    """
    if (not vcf_fn) or (not SnpPos):
        df["genotype"] = "3"
        df["Ref"] = ""
        df["Alt"] = ""
        df["ID"] = ""
        return df
    vcf_reader = pysam.VariantFile(vcf_fn)
    SnpChr, SnpCoord = SnpPos.split(":")
    SnpList = []
    hets = []
    homo_refs = []
    homo_alts = []
    SnpList = get_vcf_record(vcf_reader, SnpPos, SnpName)
    if len(SnpList) < 1:
        logging.warning("No SNP found at SnpPos")
        Ref, Alt, snpID = ["NA", "NA", "NA"]
    elif len(SnpList) > 1:
        logging.warning(
            "More than one SNP found at SnpPos. You should try to pass the argument --SnpName"
        )
        Ref, Alt, snpID = ["NA", "Ambiguous", "NA"]
    else:
        record = SnpList[0]
        Ref = record.ref
        Alt = str(record.alts[0])
        snpID = str(record.id)
        # Get sample genotypes
        hets = []
        homo_refs = []
        homo_alts = []
        for sample in record.samples:
            gt = record.samples[sample]['GT']
            # GT is a tuple, e.g. (0, 0), (0, 1), (1, 1)
            if gt is None or None in gt:
                continue
            if gt[0] == gt[1]:
                if gt[0] == 0:
                    homo_refs.append(sample)
                elif gt[0] == 1:
                    homo_alts.append(sample)
            else:
                hets.append(sample)
    df["genotype"] = df.apply(
        lambda row: label_genotype(row, homo_alts, hets, homo_refs), axis=1
    )
    df.fillna({"genotype": "3"}, inplace=True)
    df["Ref"] = Ref
    df["Alt"] = Alt
    df["ID"] = snpID
    return df


def AssignColorsFromFactorColumn(
    df,
    FactorColumnName,
    OutputColumnName="color",
    ColorPallette="RdPrBlu",
    SpecialGenotypeCase=False,
):
    # Special case for genotype column which is encoded as 0,1,2,3
    if SpecialGenotypeCase:
        ColorPallette = getattr(colorbrewer, ColorPallette)[4]
        df[OutputColumnName] = [
            rgb_to_hex(ColorPallette[int(i)]) for i in df[FactorColumnName]
        ]
    else:
        Classes = pd.factorize(df[FactorColumnName])[0]
        ColorbrewerNumClasses = max([len(set(Classes)), 3])
        ColorPallette = getattr(colorbrewer, ColorPallette)[ColorbrewerNumClasses]
        df[OutputColumnName] = [rgb_to_hex(ColorPallette[i]) for i in Classes]
    return df


def FillColorsIfNoneProvided(
    df,
    FactorColumnName,
    GroupColorColumnName="Group_color",
    OutputColumnName="color",
    ColorPallette="Dark2",
):
    df[GroupColorColumnName] = df[GroupColorColumnName].replace(
        r"^\s*$", numpy.nan, regex=True
    )
    if df[GroupColorColumnName].isnull().values.all():
        return AssignColorsFromFactorColumn(
            df=df,
            FactorColumnName=FactorColumnName,
            OutputColumnName=OutputColumnName,
            ColorPallette=ColorPallette,
        )
    else:
        df[OutputColumnName] = df[GroupColorColumnName].fillna(value="#bdbdbd")
        return df


def GetNormalizationFactorBasedOnRegionList(bigwig_fn, regionlist):
    """
    regionlist is a list of lists. Each list entries is formatted like this:
    ["chr1", 1, 10].
    """
    try:
        bw = pyBigWig.open(bigwig_fn)
        TotalSignal = []
        for (RegionChr, RegionStart, RegionStop) in regionlist:
            EntrySignal = sum(bw.values(RegionChr, RegionStart, RegionStop, numpy=True))
            TotalSignal.append(EntrySignal)
        return sum(TotalSignal) / 1e3
    except:
        return None


def NormalizeAverageAndWriteOutBigwigs(
    df,
    Region,
    dryrun=False,
    Normalization="None",
    BigwigFilesInColumnName="BigwigFilepath",
    BigwigFilesOutColumnName="bw_out",
    BigwigFilesOutPerSampleColumnName=None,
    NormalizationRegions=None,
    method="mean",
):
    """
    Given input df with column BigwigFilesInColumnName (each entry being a
    *list* of bigwig filepaths), for each row this function will write out a
    bigwig file (determined by BigwigFilesOutColumnName) where the output
    bigwig is the normalized mean value at each position. Only regions in
    Region are written out.  If Normalization=='None', no normalization occurs.
    If it is 'WholeGenome' or a filepath of bedregions, each bigwig will be
    normalized to the signal in the whole input file or the regions in
    bedregions, respectively. Also, returns the original dataframe with some
    useful columns added, like the max value after averaging (useful for
    determining plotting axis limits), and number of samples. If dryrun==T,
    then don't write out any new bigwig files, just return the dataframe with
    columns added. If BigwigFilesOutPerSampleColumnName is specified, it should
    be the name of a column in df where each row is a *list* of equal length as
    the corresponding entry BigwigFilesInColumnName, but where each entry in
    the list specifies the filepath to write out a normalized bigwig for each
    bigwig that is read in. This is required for plotting normalized coverage
    for each individual sample in addition to the default output which is just
    the group mean normalized coverage. Method can be either 'mean' or 'median'.
    """
    RegionChr, RegionCoords = Region.split(":")
    RegionStart, RegionStop = [str_to_int(i) for i in RegionCoords.split("-")]
    df["MaxAveragedValue"] = float(0)
    df["NumberSamplesAggregated"] = 0
    df["MaxPerIndValue"] = float(0)
    df["GroupsSummarisedBigwigOut"] = ""
    for i, row in df.iterrows():
        Array = []
        for j, bigwig in enumerate(row[BigwigFilesInColumnName]):
            try:
                if Normalization == "WholeGenome":
                    NormFactor = GetNormalizationFactorWholeGenome(bigwig)
                elif args.Normalization == "None":
                    NormFactor = 1
                elif args.Normalization == "PlotWindowCoverage":
                    NormFactor = GetNormalizationFactorBasedOnRegionList(bigwig, [[RegionChr, RegionStart, RegionStop]])
                else:
                    NormFactor = GetNormalizationFactorBasedOnRegionList(
                        bigwig, NormalizationRegions
                    )
                logging.debug(f"opening {bigwig}")
                bw = pyBigWig.open(bigwig)
                NormalizedValues = (
                    bw.values(RegionChr, RegionStart, RegionStop, numpy=True) / NormFactor
                )
                Array.append(NormalizedValues)
                if (not dryrun) and BigwigFilesOutPerSampleColumnName:
                    try:
                        bwOut = pyBigWig.open(
                            row[BigwigFilesOutPerSampleColumnName][j], "w"
                        )
                        bwOut.addHeader(list(bw.chroms().items()))
                        bwOut.addEntries(
                            RegionChr, RegionStart, values=NormalizedValues, span=1, step=1
                        )
                        bwOut.close()
                    except Exception as e:
                        logging.warning(f"Failed to write normalized bigwig for {bigwig}: {e}")
            except RuntimeError as e:
                logging.warning(f"Skipping bigwig {bigwig} due to error: {e}")
                continue
            except Exception as e:
                logging.warning(f"Skipping bigwig {bigwig} due to unexpected error: {e}")
                continue
        if method == "mean":
            ArrayAveraged = numpy.mean(numpy.array(Array), axis=0)
        elif method == "median":
            ArrayAveraged = numpy.median(numpy.array(Array), axis=0)
        if not dryrun:
            # write out averaged bigwig
            bwOut = pyBigWig.open(row[BigwigFilesOutColumnName], "w")
            bwOut.addHeader(list(bw.chroms().items()))
            bwOut.addEntries(
                RegionChr, RegionStart, values=ArrayAveraged, span=1, step=1
            )
            bwOut.close()
            # Save the max signal in window for to set ylimits in tracks file
        df.at[i, "MaxAveragedValue"] = max(ArrayAveraged)
        df.at[i, "NumberSamplesAggregated"] = len(Array)
        df.at[i, "MaxPerIndValue"] = numpy.max(numpy.array(Array))
        df.at[i, "GroupsSummarisedBigwigOut"] = row[BigwigFilesOutColumnName]
    return df


def GetNormalizationFactorWholeGenome(bigwig_fn):
    try:
        bw = pyBigWig.open(bigwig_fn)
        return bw.header()["sumData"] / 1e9
    except:
        return None


def WriteOutSNPBed(SnpPos, fout):
    with open(fout, "w") as f:
        if SnpPos:
            SnpChr, SnpCoord = SnpPos.split(":")
            _ = f.write("{}\t{}\t{}".format(SnpChr, str_to_int(SnpCoord) - 1, SnpCoord))


def PysamTabixToPandas_df(fn_in, Chrom=None, start=None, end=None, **kwargs):
    """
    Reading in juncs with pysam.Tabix allows memory efficient access to Region.
    To read pysam.TabixFile object into pandas, first write to StringIO object
    """
    RowsInRegion = StringIO()
    tbx = pysam.TabixFile(fn_in)
    if tbx.header:
        RowsInRegion.write(tbx.header[0] + "\n")
    for row in tbx.fetch(Chrom, start, end):
        RowsInRegion.write(row + "\n")
    RowsInRegion.seek(0)
    return pd.read_csv(RowsInRegion, sep="\t", **kwargs)


def NormalizeAverageAndWriteOutLinks(
    df,
    Region,
    FilterForJuncsBed=None,
    dryrun=False,
    BedgzFilesInColumnName="BedgzFilepath",
    SampleIDColumnName="SampleID",
    StrandColumnName="Strand",
    LinksFilesOutColumnName="links_out",
    OutputStringFormat="{0:.3g}",
    SwapStrandFilters=False,
    MinPSI=0.0,
):
    """
    Given input df with column BedgzFilesInColumn,for each row this function
    will write out a pyGenometracks links file (determined by
    LinksFilesOutColumnName) where the output links is the normalized mean
    value for each sample in SampleIDColumnName. Only junctions in Region are
    written out. Additionally, if provided, only juncs that match the chrom,
    start, stop in FilterForJuncsBed file are written out. The
    OutputStringFormat parameter can be used to help round decimals if you
    don't want sashimi plot PSI labels with too many digits. if there is a +/-
    in StrandColumnName column, only juncs on the +/- strand will be written to
    strand-specific links files. Also, returns the original dataframe with some
    useful columns added, like the max and min output value after averaging
    (useful for determining sashimi line thinkness and drawing parameters),  If
    dryrun==T, then don't write out any new bigwig files, just return the
    dataframe with columns added. If SwapStrandFilters==True, then when a row
    is marked as + strand, then only - strand juncs will be written out. This
    is useful is the biologically true strand of the bigwig and the strand
    column in the Bedgz junc files do not match. Unlike the
    NormalizeAverageAndWriteOutBigwigs function, there is no functionality for
    writing out individual links files, since I didn't think that would ever
    really be necessary. Only juncs with mean PSI > MinPSI will be written out.
    """
    df["MaxAveragedPSI"] = float(0)
    df["ContainsBedgzFile"] = "0"
    df["ContainsNonEmptyBedgzFile"] = "0"
    RegionChr, RegionCoords = Region.split(":")
    RegionStart, RegionStop = [str_to_int(i) for i in RegionCoords.split("-")]
    if SwapStrandFilters:
        PlusStrand, MinusStrand = ("-", "+")
    else:
        PlusStrand, MinusStrand = ("+", "-")
    if FilterForJuncsBed:
        try:
            JuncsToFilter_set = set(
                [
                    tuple(x)
                    for x in PysamTabixToPandas_df(
                        FilterForJuncsBed,
                        RegionChr,
                        RegionStart,
                        RegionStop,
                        usecols=[0, 1, 2],
                    ).values
                ]
            )
        except OSError:
            try:
                logging.warning(
                    f"Consider creating converting {FilterForJuncsBed} to a sorted bgziped file with tabix index for efficient access to juncs in region"
                )
                JuncsToFilter_set = set(
                    [
                        tuple(x)
                        for x in pd.read_csv(
                            FilterForJuncsBed, sep="\t", usecols=[0, 1, 2]
                        ).values
                    ]
                )
            except pd.errors.EmptyDataError:
                logging.warning(
                    "Warning, empty junc filter file... No juncs included for sashimi plot"
                )
                JuncsToFilter_set = set()
    for i, row in df.iterrows():
        if row["BedgzFilepath"]:
            df.loc[i, "ContainsBedgzFile"] = "1"
            bed_df = PysamTabixToPandas_df(
                row[BedgzFilesInColumnName], RegionChr, RegionStart, RegionStop
            )
            bed_df = bed_df.rename(
                columns={
                    bed_df.columns[0]: "#Chr",
                    bed_df.columns[1]: "start",
                    bed_df.columns[2]: "end",
                    bed_df.columns[3]: "pid",
                    bed_df.columns[5]: "Strand",
                },
                inplace=False,
            )
            bed_df["junc"] = [tuple(x) for x in bed_df.iloc[:, 0:3].values]
            # First filter by strand:
            if row["Strand"] == PlusStrand:
                bed_df = bed_df.loc[bed_df["Strand"] == "+"]
            elif row["Strand"] == MinusStrand:
                bed_df = bed_df.loc[bed_df["Strand"] == "-"]
            # Filter for juncs in FilterForJuncsBed if provided
            if FilterForJuncsBed:
                bed_df = bed_df[bed_df["junc"].isin(JuncsToFilter_set)]
            df_out = bed_df[
                bed_df.columns.intersection(
                    [str(x) for x in row[SampleIDColumnName]] + ["#Chr", "start", "end", "pid"]
                )
            ]
            df_out["psi"] = df_out.drop(["#Chr", "start", "end", "pid"], axis=1).mean(
                axis=1, skipna=True
            )
            df_out = df_out[df_out["psi"] >= float(MinPSI)]
            df_out["psi"] = df_out["psi"].apply(lambda x: OutputStringFormat.format(x))
            if len(df_out) > 0:
                df.loc[i, "ContainsNonEmptyBedgzFile"] = "1"
                df.loc[i, "MaxAveragedPSI"] = max(df_out["psi"])
            if not dryrun:
                df_out.loc[
                    :, ["#Chr", "start", "end", "#Chr", "start", "end", "psi"]
                ].to_csv(
                    row[LinksFilesOutColumnName], sep="\t", index=False, header=False
                )
    return df


def GetDefaultTemplate(DF, ColumnFactors="genotype"):
    """
    Inspect the DF and return a reasonable template from the template dir
    """
    try:
        if sys.ps1:
            DefaultTemplateDir = "./tracks_templates/"
    except:
        DefaultTemplateDir = (
            os.path.abspath(os.path.join(os.path.dirname(__file__)))
            + "/"
            + "tracks_templates/"
        )
    # if DF contains only single genotype value, then plot as
    if len(pd.factorize(DF[ColumnFactors])[1]) > 1:
        return DefaultTemplateDir + "GeneralPurposeColoredByGenotype.ini"
    else:
        return DefaultTemplateDir + "GeneralPurposeColoredByGenotype.ini"


def GetGenotypeLabel(RefString, AltString, GenotypeInt, LimitLength=5):
    GenotypeInt = int(GenotypeInt)
    if GenotypeInt == 0:
        return f"{RefString[:LimitLength]}/{RefString[:LimitLength]}"
    elif GenotypeInt == 1:
        return f"{RefString[:LimitLength]}/{AltString[:LimitLength]}"
    elif GenotypeInt == 2:
        return f"{AltString[:LimitLength]}/{AltString[:LimitLength]}"
    else:
        return ""


def get_default_output_prefix():
    tmpdir = os.environ.get("TMPDIR", "/tmp")
    return os.path.join(tmpdir, "genometracks_tmp.")

def parse_color(color_str):
    """Parse color from bed12 itemRgb field - can be RGB or hex"""
    if pd.isna(color_str) or color_str == '0' or color_str == '.':
        return '#000000'  # Default black
    
    color_str = str(color_str).strip()
    
    # Check if it's already hex
    if color_str.startswith('#'):
        return color_str
    
    # Check if it's RGB format (comma-separated)
    if ',' in color_str:
        try:
            r, g, b = [int(x.strip()) for x in color_str.split(',')]
            return f'#{r:02x}{g:02x}{b:02x}'
        except:
            return '#000000'
    
    # Try to convert single integer to RGB (assuming it's a packed RGB value)
    try:
        rgb_int = int(color_str)
        r = (rgb_int >> 16) & 0xFF
        g = (rgb_int >> 8) & 0xFF
        b = rgb_int & 0xFF
        return f'#{r:02x}{g:02x}{b:02x}'
    except:
        return '#000000'


def number_exons_by_strand(group):
    """Number exons for each transcript group based on strand orientation"""
    if group['strand'].iloc[0] == '-':
        # For minus strand, sort by start descending (3' to 5')
        group = group.sort_values('start', ascending=False)
    else:
        # For plus strand, sort by start ascending (5' to 3')
        group = group.sort_values('start', ascending=True)
    
    group['exon_number'] = range(1, len(group) + 1)
    return group


def ExtractCondensedExonCoverage(
    df,
    bed12_file,
    flanking_bp=25,
    output_prefix="",
    Region=None
):
    """
    Extract bigwig values over exon regions (plus flanking) defined in a bed12 file.
    Creates a condensed coverage matrix saved as tsv.gz for plotting.
    """
    logging.info(f"Extracting condensed exon coverage from {bed12_file}")
    
    try:
        # Get genome info from the first available bigwig
        genome_info = None
        for idx, row in df.iterrows():
            bigwig_path = row.get('bw_out', row.get('GroupsSummarisedBigwigOut', ''))
            if bigwig_path and os.path.exists(bigwig_path):
                try:
                    bw = pyBigWig.open(bigwig_path)
                    genome_info = bw.chroms()
                    bw.close()
                    logging.debug(f"Got genome info from {bigwig_path}")
                    break
                except Exception as e:
                    logging.warning(f"Could not read genome info from {bigwig_path}: {e}")
                    continue
        
        if not genome_info:
            logging.error("Could not extract genome information from any bigwig files")
            return None, None

        # Read bed12 file into pandas
        bed12_df = pd.read_csv(
            bed12_file, 
            sep='\t', 
            header=None,
            usecols=range(12),  # Ensure we only read 12 columns
            names=['chrom', 'start', 'end', 'name', 'score', 'strand', 
                   'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 
                   'blockSizes', 'blockStarts']
        )
        
        if bed12_df.empty:
            logging.warning("No bed12 entries found intersecting the specified region")
            return None, None
        
        # Create BedTool from dataframe
        bed12_bt = pybedtools.BedTool.from_dataframe(bed12_df)
        
        # Filter to region
        region_chr, region_coords = Region.split(":")
        region_start, region_stop = [str_to_int(i) for i in region_coords.split("-")]
        region_bt = pybedtools.BedTool(f"{region_chr}\t{region_start}\t{region_stop}", from_string=True)
        bed12_bt = bed12_bt.intersect(region_bt, wa=True)
        
        # Convert bed12 to bed6 (expands blocks to individual exons)
        bed6_bt = bed12_bt.bed12tobed6()
        
        # Convert bed6 BedTool to pandas DataFrame for easier manipulation
        bed6_data = []
        for exon in bed6_bt:
            bed6_data.append({
                'chrom': exon.chrom,
                'start': int(exon.start),
                'end': int(exon.end),
                'name': exon.name,
                'score': exon.score,
                'strand': exon.strand
            })
        
        if not bed6_data:
            logging.warning("No exons found after bed12 to bed6 conversion")
            return None, None
        
        bed6_df = pd.DataFrame(bed6_data)
        
        # Number exons for each transcript using the module-level function
        bed6_df = bed6_df.groupby('name').apply(number_exons_by_strand).reset_index(drop=True)
        
        # Filter to region again after conversion
        bed6_df = bed6_df[
            (bed6_df['chrom'] == region_chr) & 
            (bed6_df['start'] < region_stop) & 
            (bed6_df['end'] > region_start)
        ]
        
        # Check if any exons remain after filtering
        if bed6_df.empty:
            logging.warning("No exons found intersecting the specified region")
            return None, None
        
        # Convert back to BedTool for further processing
        bed6_bt = pybedtools.BedTool.from_dataframe(bed6_df[['chrom', 'start', 'end', 'name', 'score', 'strand']])
        
        # Create temporary genome file
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.genome', delete=False) as genome_file:
            for chrom, size in genome_info.items():
                genome_file.write(f"{chrom}\t{size}\n")
            genome_file_path = genome_file.name
        
        # Add flanking regions with slop using the temporary genome file
        bed6_flanked_bt = bed6_bt.slop(b=flanking_bp, g=genome_file_path)
        
        # Clean up temporary genome file
        os.unlink(genome_file_path)
        
        # Merge overlapping intervals
        merged_bt = bed6_flanked_bt.sort().merge()
        
        # Create intervals dataframe for plotting exons
        intervals_data = []
        
        # For each merged interval, find which original exons intersect
        for merged_interval in merged_bt:
            merged_chrom = merged_interval.chrom
            merged_start = int(merged_interval.start)
            merged_end = int(merged_interval.end)
            merged_name = f"{merged_chrom}:{merged_start}-{merged_end}"
            
            # Find intersecting original bed6 exons
            merged_interval_bt = pybedtools.BedTool(f"{merged_chrom}\t{merged_start}\t{merged_end}", from_string=True)
            intersecting_exons = bed6_bt.intersect(merged_interval_bt, wa=True)

            for exon in intersecting_exons:
                # Get exon number from our numbered dataframe
                exon_info = bed6_df[
                    (bed6_df['chrom'] == exon.chrom) & 
                    (bed6_df['start'] == int(exon.start)) & 
                    (bed6_df['end'] == int(exon.end)) & 
                    (bed6_df['name'] == exon.name)
                ]
                
                exon_number = exon_info['exon_number'].iloc[0] if len(exon_info) > 0 else 0
                
                # Get color from original bed12 entry using the module-level function
                bed12_info = bed12_df[bed12_df['name'] == exon.name]
                exon_color = '#000000'  # default
                if len(bed12_info) > 0:
                    exon_color = parse_color(bed12_info['itemRgb'].iloc[0])
                
                intervals_data.append({
                    'interval_name': merged_name,
                    'interval_start': merged_start,
                    'interval_end': merged_end,
                    'exon_chrom': exon.chrom,
                    'exon_start': int(exon.start),
                    'exon_end': int(exon.end),
                    'exon_name': exon.name if hasattr(exon, 'name') and exon.name else '',
                    'exon_strand': exon.strand if len(exon.fields) >= 6 else '.',
                    'exon_number': exon_number,
                    'exon_midpoint': int((int(exon.start) + int(exon.end)) / 2),
                    'exon_color': exon_color
                })
        
        # Create intervals dataframe and save
        if intervals_data:
            intervals_df = pd.DataFrame(intervals_data)
            intervals_df = intervals_df.drop_duplicates()  # Remove any duplicates
            intervals_df = intervals_df.sort_values(['interval_start', 'exon_start'])
            
            intervals_file = f"{output_prefix}condensed_exon_intervals.tsv.gz"
            with gzip.open(intervals_file, 'wt') as f:
                intervals_df.to_csv(f, sep='\t', index=False)
            
            logging.info(f"Exon intervals written to: {intervals_file}")
        else:
            intervals_file = None
        logging.info(f"Processed {len(bed12_df)} bed12 entries into {len(merged_bt)} merged intervals")
        
    except Exception as e:
        logging.error(f"Error processing bed12 file {bed12_file}: {e}")
        return None, None
        
    # Initialize data collection
    coverage_data = []
    
    # Process each bigwig file
    for idx, row in df.iterrows():
        bigwig_path = row.get('bw_out', row.get('GroupsSummarisedBigwigOut', ''))
        if not bigwig_path or not os.path.exists(bigwig_path):
            logging.warning(f"Bigwig file not found: {bigwig_path}")
            continue
            
        logging.debug(f"Processing bigwig: {bigwig_path}")
        
        try:
            bw = pyBigWig.open(bigwig_path)
            
            for interval in merged_bt:
                chrom = interval.chrom
                start = int(interval.start)
                end = int(interval.end)
                interval_name = f"{chrom}:{start}-{end}"
                
                # Extract values
                try:
                    values = bw.values(chrom, start, end, numpy=True)
                    if values is not None and len(values) > 0:
                        # Create position array
                        positions = numpy.arange(start, end)
                        
                        # Store data for each position
                        for pos, val in zip(positions, values):
                            if not numpy.isnan(val):  # Skip NaN values
                                # Start with the basic coverage info
                                data_row = {
                                    'chrom': chrom,
                                    'position': pos,
                                    'interval_name': interval_name,
                                    'interval_start': start,
                                    'interval_end': end,
                                    'position_in_interval': pos - start,
                                    'coverage': val,
                                    'bigwig_file': os.path.basename(bigwig_path),
                                }
                                
                                # Add ALL columns from the DF row (except list columns which need special handling)
                                for col in df.columns:
                                    if col not in ['BigwigFilepath', 'bw_out_PerInd', 'SampleID']:  # Skip list columns
                                        data_row[col] = row[col]
                                    elif col == 'SampleID':
                                        # Convert SampleID list to comma-separated string
                                        data_row['sample_ids'] = ','.join(row['SampleID']) if isinstance(row['SampleID'], list) else str(row['SampleID'])
                                
                                coverage_data.append(data_row)
                                
                except Exception as e:
                    logging.warning(f"Error extracting values from {chrom}:{start}-{end}: {e}")
                    continue
                    
            bw.close()
            
        except Exception as e:
            logging.warning(f"Error processing bigwig {bigwig_path}: {e}")
            continue
    
    # Convert to DataFrame and save
    coverage_file = None
    if coverage_data:
        coverage_df = pd.DataFrame(coverage_data)
        
        # Sort by genomic position
        coverage_df = coverage_df.sort_values(['chrom', 'position', 'Group_label', 'genotype'])
        
        # Output file
        coverage_file = f"{output_prefix}condensed_exon_coverage.tsv.gz"
        
        # Write compressed TSV
        with gzip.open(coverage_file, 'wt') as f:
            coverage_df.to_csv(f, sep='\t', index=False)
            
        logging.info(f"Condensed exon coverage written to: {coverage_file}")
        logging.info(f"Coverage data shape: {coverage_df.shape}")
        logging.info(f"Intervals covered: {len(coverage_df['interval_name'].unique())}")
        logging.info(f"Columns in output: {list(coverage_df.columns)}")
    else:
        logging.warning("No coverage data extracted")
    
    return coverage_file, intervals_file

def plot_condensed_exon_coverage(coverage_file, intervals_file, output_pdf, reverse_xaxis=False, figsize=(16, 10), dpi=300):
    """
    Create a condensed exon coverage plot similar to the R ggplot2 version.
    
    Parameters:
    -----------
    coverage_file : str
        Path to condensed_exon_coverage.tsv.gz file
    intervals_file : str
        Path to condensed_exon_intervals.tsv.gz file
    output_pdf : str
        Output PDF file path
    reverse_xaxis : bool
        Whether to reverse the x-axis and facet column order.
    figsize : tuple
        Figure size (width, height) in inches
    dpi : int
        Resolution for the output
    """
    from matplotlib.backends.backend_pdf import PdfPages
    
    # Clear any matplotlib state
    plt.close('all')
    
    # Read the data - force fresh read
    logging.info(f"Reading coverage data from {coverage_file}")
    dat = pd.read_csv(coverage_file, sep='\t', compression='gzip')
    dat = dat.copy()  # Force a fresh copy
    
    logging.info(f"Reading intervals data from {intervals_file}")
    exons = pd.read_csv(intervals_file, sep='\t', compression='gzip')
    exons = exons.copy()  # Force a fresh copy
    
    # Transform coverage values: make negative for minus strand
    dat['coverage_transformed'] = dat['coverage'].copy()
    dat.loc[dat['Strand'] == '-', 'coverage_transformed'] = -dat.loc[dat['Strand'] == '-', 'coverage']
    
    # Debug: Print what we actually read
    logging.debug(f"Coverage data columns: {list(dat.columns)}")
    logging.debug(f"Unique Supergroups in data: {sorted(dat['Supergroup'].unique())}")
    logging.debug(f"Unique Strands in data: {sorted(dat['Strand'].unique())}")
    logging.debug(f"Data shape: {dat.shape}")
    
    # Sort FullLabel by PlotOrder for correct plotting order
    plot_order = dat.drop_duplicates(['FullLabel', 'PlotOrder']).sort_values('PlotOrder')
    fulllabel_order = plot_order['FullLabel'].tolist()
    dat['FullLabel'] = pd.Categorical(dat['FullLabel'], categories=fulllabel_order, ordered=True)
    
    # Create color mapping
    colors_dict = dat.drop_duplicates(['FullLabel', 'color']).set_index('FullLabel')['color'].to_dict()
    
    # Get unique combinations for faceting - use Supergroup instead of SupergroupLabel
    supergroups = sorted(dat['Supergroup'].unique())
    
    # Debug: Confirm what we're using
    logging.info(f"Using supergroups: {supergroups}")
    
    # Sort intervals based on reverse_xaxis flag
    if reverse_xaxis:
        intervals = exons.sort_values('interval_start', ascending=False)['interval_name'].unique()
    else:
        intervals = sorted(dat['interval_name'].unique())
    
    logging.info(f"Using intervals: {list(intervals)}")
    
    # Calculate widths for each interval based on genomic coordinates
    interval_widths = []
    interval_info = {}
    for interval in intervals:
        interval_data = dat[dat['interval_name'] == interval]
        if len(interval_data) > 0:
            x_min = interval_data['position'].min()
            x_max = interval_data['position'].max()
            width = x_max - x_min
            interval_widths.append(width)
            interval_info[interval] = {'start': x_min, 'end': x_max, 'width': width}
        else:
            interval_widths.append(1)  # Default width for empty intervals
            interval_info[interval] = {'start': 0, 'end': 1, 'width': 1}
    
    # Normalize widths to sum to 1 for GridSpec width_ratios
    total_width = sum(interval_widths)
    width_ratios = [w / total_width for w in interval_widths]
    
    logging.info(f"Interval widths: {dict(zip(intervals, interval_widths))}")
    logging.info(f"Width ratios: {dict(zip(intervals, width_ratios))}")
    
    # Calculate y-axis limits for each supergroup (row)
    supergroup_ylims = {}
    for supergroup in supergroups:
        supergroup_data = dat[dat['Supergroup'] == supergroup]
        if len(supergroup_data) > 0:
            y_min = supergroup_data['coverage_transformed'].min()
            y_max = supergroup_data['coverage_transformed'].max()
            # Add some padding
            y_range = y_max - y_min
            y_padding = y_range * 0.1 if y_range > 0 else 1
            supergroup_ylims[supergroup] = (y_min - y_padding, y_max + y_padding)
        else:
            supergroup_ylims[supergroup] = (0, 1)
    
    # Create the plot with even more space for legend
    fig = plt.figure(figsize=(figsize[0] + 4, figsize[1]), dpi=dpi)
    
    # Calculate subplot layout
    n_supergroups = len(supergroups)
    n_intervals = len(intervals)
    
    logging.info(f"Creating plot with {n_supergroups} supergroups and {n_intervals} intervals")
    
    # Create subplots with width_ratios for proportional spacing
    gs = fig.add_gridspec(n_supergroups, n_intervals, 
                         width_ratios=width_ratios,  # This is the key addition!
                         hspace=0.2, wspace=0.06,
                         left=0.15, right=0.7, top=0.95, bottom=0.15)
    
    # Track which axes have data for legend creation
    axes_with_data = []
    
    # Plot for each combination of supergroup and interval
    for i, supergroup in enumerate(supergroups):
        logging.debug(f"Processing supergroup {i}: {supergroup}")
        for j, interval in enumerate(intervals):
            ax = fig.add_subplot(gs[i, j])
            
            # Filter data for this facet - use Supergroup instead of SupergroupLabel
            facet_data = dat[(dat['Supergroup'] == supergroup) & 
                           (dat['interval_name'] == interval)]
            
            logging.debug(f"Facet [{i},{j}] ({supergroup}, {interval}): {len(facet_data)} data points")
            
            if len(facet_data) > 0:
                # Plot lines for each FullLabel using transformed coverage values
                for full_label in facet_data['FullLabel'].unique():
                    label_data = facet_data[facet_data['FullLabel'] == full_label]
                    if len(label_data) > 0:
                        ax.plot(label_data['position'], label_data['coverage_transformed'], 
                               color=colors_dict.get(full_label, '#000000'),
                               label=full_label, linewidth=1.5)
                
                # Set x-axis limits with some padding
                x_min, x_max = facet_data['position'].min(), facet_data['position'].max()
                x_padding = (x_max - x_min) * 0.02
                ax.set_xlim(x_min - x_padding, x_max + x_padding)
                
                # Set y-axis limits to match all facets in this row (supergroup)
                y_min_row, y_max_row = supergroup_ylims[supergroup]
                ax.set_ylim(y_min_row, y_max_row)
                
                # Add exon segments ONLY to bottom row
                if i == n_supergroups - 1:  # Only bottom row
                    facet_exons = exons[exons['interval_name'] == interval]
                    if len(facet_exons) > 0:
                        # Position exon segments at the bottom of the plot, but within the visible area
                        y_min_row, y_max_row = supergroup_ylims[supergroup]
                        y_range = y_max_row - y_min_row
                        y_bottom = y_min_row + (y_range * 0.02)  # Position just above the bottom edge
                        
                        for _, exon_row in facet_exons.iterrows():
                            # Use the exon's color if available, otherwise default to black
                            exon_color = exon_row.get('exon_color', '#000000')
                            ax.hlines(y=y_bottom, xmin=exon_row['exon_start'], xmax=exon_row['exon_end'],
                                     colors=exon_color, linewidth=3)
                
                # Add horizontal line at y=0 to separate positive and negative values
                ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.5, alpha=0.7)
                
                # Formatting - equivalent to theme_classic()
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.spines['left'].set_linewidth(0.5)
                ax.spines['bottom'].set_linewidth(0.5)
                
                # Reverse x-axis if the flag is set
                if reverse_xaxis:
                    ax.invert_xaxis()
                
                # Only show x-axis labels and ticks on BOTTOM ROW
                if i == n_supergroups - 1:  # Bottom row
                    # Use exon-based ticks and labels
                    facet_exons = exons[exons['interval_name'] == interval]
                    if len(facet_exons) > 0:
                        # Get exon midpoints and numbers for ticks
                        exon_ticks = facet_exons['exon_midpoint'].tolist()
                        exon_labels = [f"E{int(num)}" for num in facet_exons['exon_number'].tolist()]
                        
                        ax.set_xticks(exon_ticks)
                        ax.set_xticklabels(exon_labels, rotation=0, ha='center')
                    else:
                        # Fallback to coordinate-based ticks if no exon info
                        x_ticks = numpy.linspace(x_min, x_max, 3)
                        ax.set_xticks(x_ticks)
                        ax.set_xticklabels([f'{int(tick):,}' for tick in x_ticks], 
                                          rotation=45, ha='right')
                else:  # Not bottom row - remove x-axis labels
                    ax.set_xticks([])
                    ax.set_xticklabels([])
                    ax.spines['bottom'].set_visible(False)
                
                # Only show y-axis labels and ticks on LEFTMOST COLUMN
                if j == 0:  # Leftmost column
                    # Keep default y-axis labels
                    pass
                else:  # Not leftmost column - remove y-axis labels
                    ax.set_yticks([])
                    ax.set_yticklabels([])
                    ax.spines['left'].set_visible(False)
                
                # Remove axis labels and titles
                ax.set_xlabel('')
                ax.set_ylabel('')
                ax.set_title('')
                
                # Store axis for legend if it has data
                if ax.get_lines():
                    axes_with_data.append(ax)
                
            else:
                # Empty facet - set same y-limits as other facets in this row and hide it
                y_min_row, y_max_row = supergroup_ylims[supergroup]
                ax.set_ylim(y_min_row, y_max_row)
                ax.set_visible(False)
            
            # Add Supergroup text on the leftmost column
            if j == 0:  # Leftmost column
                # Calculate the position for the text (middle of the subplot)
                bbox = ax.get_position()
                y_center = bbox.y0 + bbox.height / 2
                
                # Add text closer to the axis (between axis and Coverage label)
                logging.debug(f"Adding supergroup label '{supergroup}' at position {y_center}")
                fig.text(0.10, y_center, supergroup, 
                        rotation=90, va='center', ha='center', 
                        fontsize=10, weight='bold')
    
    # Add overall labels
    fig.text(0.5, 0.05, 'Exon Number', ha='center', fontsize=12, weight='bold')
    fig.text(0.04, 0.5, 'Coverage', va='center', rotation='vertical', fontsize=12, weight='bold')
    
    # Create legend using the first subplot that has data, ordered by PlotOrder
    if axes_with_data:
        # Get all unique labels from the plot and sort by PlotOrder
        all_labels_from_plot = set()
        for ax in axes_with_data:
            handles, labels = ax.get_legend_handles_labels()
            all_labels_from_plot.update(labels)
        
        # Create ordered legend based on PlotOrder
        ordered_handles = []
        ordered_labels = []
        
        for full_label in fulllabel_order:
            if full_label in all_labels_from_plot:
                # Create a dummy line with the right color for the legend
                dummy_line = plt.Line2D([0], [0], color=colors_dict.get(full_label, '#000000'), linewidth=1.5)
                ordered_handles.append(dummy_line)
                ordered_labels.append(full_label)
        
        # Place legend in the dedicated space on the right with more room
        fig.legend(ordered_handles, ordered_labels, 
                  loc='center left', bbox_to_anchor=(0.72, 0.5),
                  frameon=False, fontsize=10)
    
    # Save to PDF
    with PdfPages(output_pdf) as pdf:
        pdf.savefig(fig, bbox_inches='tight', dpi=dpi)
        logging.info(f"Plot saved to {output_pdf}")
    
    plt.close(fig)
    return output_pdf

def parse_args(Args=None):
    p = argparse.ArgumentParser(
        conflict_handler="resolve",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--BigwigList",
        help=""" If using with the --BigwigListType KeyFile  option, a tab delimited text file with samples in the first column and a path to the input bigwig files in the second column. Alternatively, if using with the --BigwigListType GlobPattern option, use a  wildcard glob of bigwig files, and the samples will be inferred by searching the expanded filepaths for strings that match samples in the VCF. Note that using this wildcard method requires you to enclose the grob pattern in quotes, since this script needs to parse the grob pattern with python. Excluding quotes will immediately expand the grob pattern by the shell.  Example1: MySampleToBigwigKey.tsv --BigwigListType KeyFile Example2: "./Bigwigs/*/Coverage.bw" --BigwigListType GlobPattern where the * expansion may be sample names, which are automatically matched (with some wiggle room). Using a keyfile is obviously more robust, as the glob pattern requires some guessing to match files to samples. For convenient use with 1000Genomes project samples, this sripts will match samples by searching filenames for 5 digit patterns that match sample naming scheme used in 1000Genomes project. For example: a file named ./MyBigwigs/GM19137.bw will be matched to the sample HG19137 because of the 19137 substring match.  If a KeyFile is provided, it should be a 2+ column tab delimited file where the first 5 columns defined as: 1: sample ID (must match vcf sample if attempting to average by genotype) 2: bigwig filepath 3: group_label (optional) 4: ""|"+"|"-" for strand. (Do not actually enter quotes.. optional, default is "" for unstranded. this is useful for plotting stranded coverage.  Columns past these 4 will be ignored""",
        metavar='{"<GlobPattern>",<KEYFILE>}',
    )
    p.add_argument(
        "--Region",
        help="Region to output",
        metavar="<CHR>:<START>-<STOP>",
        required=True,
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
        "--BigwigListType",
        choices=["KeyFile", "GlobPattern"],
        help="Required. Define how the bigwig list positional argument is specified. The GlobBattern option requires use of --VCF option",
        required=True,
    )
    p.add_argument(
        "--Normalization",
        help="A bed file of regions to use for sample-depth normalization. The within-sample total coverage over the regions will be used to normalize each sample. If 'None' is chosen, will not perform any normalization. If 'WholeGenome' is used, will use whole genome. 'PlotWindowCoverage' normalizes tracks to total coverage specified in --Region. I haven't ever tested the <BEDFILE> option. (default: %(default)s)",
        default="None",
        metavar="{<BEDFILE>,None,WholeGenome,PlotWindowCoverage}",
    )
    p.add_argument(
        "--OutputPrefix",
        default=get_default_output_prefix(),
        help="Prefix for output files (default: a temp directory, usually $TMPDIR or /tmp)."
    )
    p.add_argument(
        "--BedfileForSashimiLinks",
        help="QTLtools style bed or bed.gz file with header for samples (sample names must match vcf) and values to calculate average PSI per genotype. (default: %(default)s)",
        default=None,
    )
    p.add_argument(
        "--GroupSettingsFile",
        help=""" Requires use of --BigwigListType KeyFile. This option specifies tab delimited file with one row per Group_label to add extra plotting settings for each group. Columns are as follows: 1: Group_label (must match a Group_label in the KeyFile) 2: Group_color (optional). Hex or rgb colors to plot potentially plot each group as defined in the output ini 3: BedgzFile (optional). Any bedgz file provided in this column will take precedent over bedgz files provided to the --BedfileForSashimiLinks option. If a group settings file is provided only groups described in that file will be plotted. Additionally, tracks will be plotted in the order they are listed in this file""",
        default=None,
    )
    p.add_argument(
        "--OutputNormalizedBigwigsPerSample",
        help="Output normalized bigwigs for each sample. This will be required if the ini template chosen plots coverage of individual samples.",
        action="store_true",
    )
    p.add_argument(
        "--TracksTemplate",
        metavar="<FILE>",
        help="A jinja template for a tracks file for pyGenomeTracks customization. An example is included. Template variables allowed are 'OutputPrefix', 'HomoRefTitle', 'HetTitle', 'HomoAltTitle', and 'YMax'. If this argument is provided, the template file will be populated with the template variables to create a tracks file that can be used for pyGenomeTracks. If this argument is not provided, will output a very basic tracks file that can be used for pyGenomeTracks",
        default=None,
    )
    p.add_argument(
        "--IGVSessionOutput",
        action="store_true",
        help="Flag to create a tracks.xml IGV session file (see https://software.broadinstitute.org/software/igv/Sessions) output to easily browse aggregated bigwigs with IGV.",
        default=False,
    )
    p.add_argument(
        "--JinjaDfOutput",
        action="store_true",
        help="Flag to create a DF.tsv file to see the dataframe object (accessible in jinja template as variable named 'DF') used to create the ini file",
        default=False,
    )
    p.add_argument(
        "--IGVSessionTemplate",
        metavar="<FILE>",
        help="A jinja template for a IGV session xml file. The template file will be populated with the DF variables to create the optional output file specified by IGVSessionOutput. If this argument is not provided, will basic template will be chosen",
        default=None,
    )
    p.add_argument(
        "--FilterJuncsByBed",
        metavar="<FILE>",
        help="An optional bedfile of junctions to filter for inclusion. If none is provided, all juncs will be included",
        default=None,
    )
    p.add_argument(
        "--NoSashimi",
        help="Flag to disable sashimi plots. Default is to include sashimi plot if junction file(s) are provided (as a column in a KeyFile, or with the --BedfileForSashimiLinks argument) and there are juncs in the plotting region. This flag disables sashimi plots regardless",
        action="store_true",
        default=False,
    )
    p.add_argument(
        "--Bed12GenesToIni",
        metavar="<FILE>",
        help="An optional bed12 genes files. If provided, will add a track for these genes at the end of the ini output file",
        default=None,
    )
    p.add_argument(
        "--GenotypeBrewerPalettes",
        metavar="<string>",
        help=f'A colorbrewer palette to use for genotypes. Options: {colorbrewer_palletes}',
        default='RdPrBlu',
    )
    p.add_argument(
        "--GenotypeCustomPalette",
        metavar="<#NNNNNN,#NNNNNN,#NNNNNN,#NNNNNN>",
        help='string of 4 commad delimited HEX color codes to use for genotypes: (1) ref/ref, (2) ref/alt, (3) alt/alt, and (4) ungenotyped. If provided, this argument takes precedence over --GenotypeBrewerPalette',
        default=None,
    )
    p.add_argument(
        "--Workdir",
        metavar="<path>",
        help="An optional path to set as workdir before executing main script function. Could be useful, for example, if bigwig file list uses relative file paths from a different directory",
        default="./",
    )
    p.add_argument(
        "--pyGenomeTracks_args",
        nargs=argparse.REMAINDER,
        help=(
            "Arguments to pass directly to pyGenomeTracks after ini file is generated. "
            "By default, '--region <Region>' will be appended unless you specify a --region argument here. "
            "Usage: --pyGenomeTracks_args <args for pyGenomeTracks>"
        ),
    )
    p.add_argument(
        "--CondensedExonViewBed",
        metavar="<FILE>",
        help="A bed12 file containing blocked entries (genes/transcripts with exons). If provided, will output a condensed coverage matrix (tsv.gz) containing bigwig values only over exonic regions plus flanking sequence for compact visualization.",
        default=None,
    )
    p.add_argument(
        "--CondensedExonView_DecreasingXAxis",
        action="store_true",
        help="If provided, reverse the x-axis and facet column order for minus strand genes in the condensed exon view plot.",
        default=False,
    )
    p.add_argument(
        "--CondensedExonView_FigSize",
        nargs=2,
        type=float,
        default=[16.0, 10.0],
        metavar=('WIDTH', 'HEIGHT'),
        help="Figure size (width, height) in inches for condensed exon view plot (default: %(default)s)",
    )
    p.add_argument(
        "--ExonFlankingBp",
        type=int,
        default=25,
        help="Number of base pairs of flanking intronic sequence to include around each exon for condensed view (default: %(default)s)",
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


def main(**kwargs):
    """
    return dataframe of samples, bigwig files, and do main work, of this script, resulting in output files including a tracks.ini file for pygenometracks.
    """
    logging.debug(kwargs)
    os.chdir(kwargs["Workdir"])
    logging.info("Reading in bigwiglist")
    if kwargs["BigwigListType"] == "KeyFile":
        logging.debug("readKeyFile")
        DF = ReadSamplesTsv(kwargs["BigwigList"])
    elif kwargs["BigwigListType"] == "GlobPattern":
        if kwargs["VCF"]:
            logging.info("Making Key File")
            DF = ConstructSamplesTsvFromGlobPattern(kwargs["BigwigList"], kwargs["VCF"])
        else:
            raise Exception("Files bigwig files as a glob pattern, must supply a VCF")
    logging.info("Adding genotype to DF")
    DF = AppendGenotypeColumnToDf(DF, vcf_fn=kwargs["VCF"], SnpPos=kwargs["SnpPos"], SnpName=kwargs["SnpName"])
    logging.info("Adding column for individual output bigwigs filepaths to DF")
    DF["bw_out_PerInd"] = (
        kwargs["OutputPrefix"]
        + DF[["SampleID", "Group_label", "genotype", "Strand"]].astype(str).agg("-".join, axis=1)
        + ".bw"
    )
    # groubpy
    DF = DF.groupby(
        [
            col
            for col in DF.columns
            if col not in ["BigwigFilepath", "SampleID", "bw_out_PerInd"]
        ],
        as_index=False,
    )[["BigwigFilepath", "bw_out_PerInd", "SampleID"]].agg(lambda x: list(x))
    logging.info("Mapping group setting to samples")
    DF = AddExtraColumnsPerGroup(
        DF,
        GroupToBedTsv=kwargs["GroupSettingsFile"],
        BedfileForAll=kwargs["BedfileForSashimiLinks"],
    )
    # Add column for group mean bigwig filepaths
    DF["bw_out"] = (
        kwargs["OutputPrefix"]
        + DF[["Group_label", "genotype", "Strand"]].agg("-".join, axis=1)
        + ".bw"
    )
    # Add column for genotype colors
    if kwargs['GenotypeCustomPalette']:
        colorbrewer.Custom = {
            4: [hex_to_rgb(c) for c in kwargs['GenotypeCustomPalette'].split(',')],
        }
        kwargs['GenotypeBrewerPalettes'] = 'Custom'
    DF = AssignColorsFromFactorColumn(
        DF, "genotype", "genotype_color", SpecialGenotypeCase=True, ColorPallette=kwargs['GenotypeBrewerPalettes']
    )
    # Add column for group color, taking from existing column if not empty
    DF = FillColorsIfNoneProvided(
        DF,
        "Group_label",
        GroupColorColumnName="Group_color",
        OutputColumnName="Group_color_final",
    )
    # write out mean bw files, and normalized individual bws if in necessary
    if kwargs["OutputNormalizedBigwigsPerSample"]:
        BigwigFilesOutPerSampleColumnName = "bw_out_PerInd"
    else:
        BigwigFilesOutPerSampleColumnName = None
    logging.info("Writing normalized and averaged bigwigs")
    DF = NormalizeAverageAndWriteOutBigwigs(
        DF,
        kwargs["Region"],
        dryrun=False,
        Normalization=kwargs["Normalization"],
        BigwigFilesOutPerSampleColumnName=BigwigFilesOutPerSampleColumnName,
    )
    # Add column for links files output and write out links files
    DF["links_out"] = (
        kwargs["OutputPrefix"]
        + "MeanPSI."
        + DF[["Group_label", "genotype", "Strand"]].agg("-".join, axis=1)
        + ".links"
    )
    DF["ContainsBedgzFile"] = "0"
    DF["ContainsNonEmptyBedgzFile"] = "0"
    if not kwargs["NoSashimi"]:
        # Check if any bed files are available for sashimi plotting
        has_bedfiles = False
        
        # Check if BedfileForSashimiLinks is provided
        if kwargs["BedfileForSashimiLinks"]:
            has_bedfiles = True
        
        # Check if any groups have bedgz files in the BedgzFilepath column
        if not has_bedfiles and 'BedgzFilepath' in DF.columns:
            # Check for non-empty bedgz file paths
            non_empty_bedgz = DF['BedgzFilepath'].fillna('').str.strip()
            if (non_empty_bedgz != '').any():
                has_bedfiles = True
        
        if has_bedfiles:
            logging.info("Writing and averaged links")
            DF = NormalizeAverageAndWriteOutLinks(
                DF,
                kwargs["Region"],
                dryrun=False,
                FilterForJuncsBed=kwargs["FilterJuncsByBed"],
                MinPSI=1,
            )
        else:
            logging.warning(
                "No bed files provided for sashimi links. Either provide --BedfileForSashimiLinks, "
                "add BedgzFilepath entries to your GroupSettingsFile, or use --NoSashimi to disable "
                "sashimi plotting entirely."
            )
            # Set default values for sashimi-related columns to prevent downstream errors
            DF["MaxAveragedPSI"] = 0.0
            DF["ContainsBedgzFile"] = "0"
            DF["ContainsNonEmptyBedgzFile"] = "0"
    else:
        logging.info("Sashimi plotting disabled by --NoSashimi flag")
        # Set default values for sashimi-related columns
        DF["MaxAveragedPSI"] = 0.0
        DF["ContainsBedgzFile"] = "0" 
        DF["ContainsNonEmptyBedgzFile"] = "0"
    # Add some more columns that could be useful for ini template
    DF["PerGroupMaxMean"] = DF.groupby(["Group_label"])["MaxAveragedValue"].transform(
        max
    )
    DF["PerGroupMaxPerInd"] = DF.groupby(["Group_label"])["MaxPerIndValue"].transform(
        max
    )
    DF["PerSupergroupMaxMean"] = DF.groupby(["Supergroup"])["MaxAveragedValue"].transform(
        max
    )
    DF["PerSupergroupMaxPerInd"] = DF.groupby(["Supergroup"])["MaxPerIndValue"].transform(
        max
    )
    DF['NumGroupsInSupergroup'] = DF.groupby('Supergroup')['Supergroup'].transform('count')
    if len(set(DF["genotype"])) == 1:
        DF["color"] = DF["Group_color_final"]
    else:
        DF["color"] = DF["genotype_color"]
    # Add some extra columns that could be useful labels in jinja ini template
    DF["genotype_label"] = DF.apply(
        lambda row: GetGenotypeLabel(row["Ref"], row["Alt"], row["genotype"]), axis=1
    )
    DF["ShortLabel"] = [
        f"{a} ({b}) {c}"
        for a, b, c in zip(
            DF["Group_label"], DF["NumberSamplesAggregated"], DF["Strand"]
        )
    ]
    DF["FullLabel"] = [
        f"{a} {b} ({c}) {d}"
        for a, b, c, d in zip(
            DF["Group_label"],
            DF["genotype_label"],
            DF["NumberSamplesAggregated"],
            DF["Strand"],
        )
    ]
    DF["SupergroupLabel"] = [
        f"{a} ({b}) {c}"
        for a, b, c in zip(
            DF["Supergroup"],
            DF["NumGroupsInSupergroup"],
            DF["Strand"],
        )
    ]
    MockAxisLabels_DF = DF.groupby(["Group_label", "Strand"], as_index=False)[
        ["genotype", "NumberSamplesAggregated"]
    ].agg(lambda x: tuple(x))
    DF = DF.merge(
        MockAxisLabels_DF,
        how="left",
        on=["Group_label", "Strand"],
        suffixes=["", "_LabelForMockAxes"],
    )
    DF["MockAxesFullLabel"] = [
        f"{a} {b} {c}"
        for a, b, c in zip(
            DF["Group_label"],
            DF["NumberSamplesAggregated_LabelForMockAxes"],
            DF["Strand"],
        )
    ]
    WriteOutSNPBed(kwargs["SnpPos"], kwargs["OutputPrefix"] + "SNP.bed")
    DF = DF.sort_values(
        by=["PlotOrder", "Supergroup", "Group_label", "Strand", "genotype"], ascending=[True, True, True, True, False]
    ).reset_index(drop=True)
    # import pdb; pdb.set_trace()
    # Get jinja2 template
    if kwargs["TracksTemplate"]:
        with open(kwargs["TracksTemplate"], "r") as fh:
            template = Template(fh.read())
    else:
        with open(GetDefaultTemplate(DF), "r") as fh:
            template = Template(fh.read())
    logging.info(
        "Writing template with dataframe with the following df.columns accessible in jinja2 template: {}".format(
            DF.columns
        )
    )
    with open(kwargs["OutputPrefix"] + "tracks.ini", "w") as f:
        _ = f.write(template.render(DF=DF, OutputPrefix=kwargs["OutputPrefix"], Bed12GenesFile=kwargs["Bed12GenesToIni"]))
    # write IGV xml
    # import pdb; pdb.set_trace()
    if kwargs["JinjaDfOutput"]:
        DF.to_csv(kwargs["OutputPrefix"] + "DF.tsv", sep="\t")
    if kwargs["IGVSessionOutput"]:
        if kwargs["IGVSessionTemplate"]:
            with open(kwargs["IGVSessionTemplate"], "r") as fh:
                template = Template(fh.read())
        else:
            with open("/project2/yangili1/bjf79/20211209_JingxinRNAseq/code/scripts/GenometracksByGenotype/tracks_templates/GeneralPurpose.xml", "r") as fh:
                template = Template(fh.read())
        logging.info(
            "Writing xml template with dataframe with the following df.columns accessible in jinja2 template: {}".format(
                DF.columns
            )
        )
        with open(kwargs["OutputPrefix"] + "tracks.xml", "w") as f:
            _ = f.write(template.render(DF=DF, Region=args.Region, OutputPrefix=kwargs["OutputPrefix"], Bed12GenesFile=kwargs["Bed12GenesToIni"]))
    if kwargs.get("pyGenomeTracks_args"):
        import subprocess
        ini_file = kwargs["OutputPrefix"] + "tracks.ini"
        pygt_args = kwargs["pyGenomeTracks_args"]
        # Check if '--region' is already in the args
        if not any(arg == "--region" or arg.startswith("--region=") for arg in pygt_args):
            pygt_args = ["--region", kwargs["Region"]] + pygt_args
        cmd = ["pyGenomeTracks", "--tracks", ini_file] + pygt_args
        print("Running:", " ".join(cmd))
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            logging.info(result.stdout)
        except subprocess.CalledProcessError as e:
            logging.error(f"pyGenomeTracks failed: {e.stderr}")
            raise

    if kwargs.get("CondensedExonViewBed"):
        logging.info("Generating condensed exon coverage matrix")
        condensed_file, intervals_file = ExtractCondensedExonCoverage(
            DF,
            kwargs["CondensedExonViewBed"],
            flanking_bp=kwargs.get("ExonFlankingBp", 25),
            output_prefix=kwargs["OutputPrefix"],
            Region=kwargs.get("Region"))
        if condensed_file:
            logging.info(f"Condensed coverage matrix saved to: {condensed_file}")
        if intervals_file:
            logging.info(f"Exon intervals saved to: {intervals_file}")
        # Generate the plot
        if condensed_file and intervals_file:
            plot_output = kwargs["OutputPrefix"] + "condensed_exon_coverage.pdf"
            try:
                figsize = tuple(kwargs.get("CondensedExonView_FigSize", [16.0, 10.0]))
                plot_condensed_exon_coverage(
                    condensed_file, 
                    intervals_file, 
                    plot_output, 
                    reverse_xaxis=kwargs.get("CondensedExonView_DecreasingXAxis", False),
                    figsize=figsize
                )
                logging.info(f"Condensed exon coverage plot saved to: {plot_output}")
            except Exception as e:
                logging.warning(f"Could not generate condensed exon coverage plot: {e}")
    return DF

if __name__ == "__main__":
    # I like to script and debug with an interactive interpreter in the same
    # working dir as this script.. If using interactive interpreter, can step
    # thru the script with args defined below
    try:
        if sys.ps1:
            # test two groups, by genotype, with strandedness
            # Args = '--VCF /project2/yangili1/bjf79/ChromatinSplicingQTLs/code/Genotypes/1KG_GRCh38_SubsetYRI/2.vcf.gz --SnpPos chr2:74482972  --Normalization None --BigwigListType KeyFile --OutputPrefix test_results/ -vvvv --BigwigList test_data/sample_list.tsv --Region chr2:74,480,713-74,505,757'.split(' ')
            # Test splicing by genotype
            ### NOTE that glob pattern should be quoted to prevent bash expansion when args parsed from command line. When args parsed from python, do not quote
            # Args = '--Normalization None --BigwigListType KeyFile --OutputPrefix test_results/ -vv --BigwigList test_data/sample_list.tsv --Region chr2:74,480,713-74,505,757'.split(' ')
            # Args = "--Workdir /project2/yangili1/bjf79/ChromatinSplicingQTLs/code --Region chr1:209745471-209812326 --SnpPos chr1:209806682 --VCF Genotypes/1KG_GRCh38/1.vcf.gz --BigwigListType KeyFile --GroupSettingsFile PlotQTLs/bwList.Groups.tsv --BigwigList PlotQTLs/bwList.tsv --OutputPrefix scratch/testplot.chunk1. --TracksTemplate /project2/yangili1/bjf79/GenometracksByGenotype/tracks_templates/GeneralPurposeColoredByGenotype.ini -q".split(" ")
            # Test bug for Luna
            Args = "--BigwigList test_data/bigwig.list.OnlyFirstExample.tsv --BigwigListType KeyFile --Region chr4:3,211,727-3,215,527 --GroupSettingsFile test_data/Example_HTT_data/GroupsFile.A.tsv --Bed12GenesToIni PremadeTracks/gencode.v26.FromGTEx.genes.bed12.gz --BedfileForSashimiLinks test_data/Example_HTT_data/SplicingTable.bed.gz --pyGenomeTracks_args -o test.pdf".split(" ")
            Args = "--BigwigList /project2/yangili1/zeqingj/files/sashimi/bwList.fixed.tsv --BigwigListType KeyFile --VCF /project2/yangili1/bjf79/ChromatinSplicingQTLs/code/Genotypes/1KG_GRCh38/Autosomes.vcf.gz --SnpPos chr7:128675737 --Region chr7:128573791-128676737 --GroupSettingsFile /project2/yangili1/zeqingj/files/sashimi/bwList.Groups.4col.tsv --NoSashimi --pyGenomeTracks_args -o test.pdf".split(" ")
            args = parse_args(Args=Args)
    except:
        args = parse_args()
    setup_logging(args.verbosity)
    DF = main(**vars(args))
    DF.to_csv(args.OutputPrefix + "DF.final.tsv", sep="\t", index=False)




