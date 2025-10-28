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
import argparse
import logging

def setup_logging(verbosity):
    """Set up logging based on verbosity level."""
    if verbosity == 1:
        level = logging.WARNING
    elif verbosity == 2:
        level = logging.INFO
    elif verbosity >= 3:
        level = logging.DEBUG
    else:
        level = logging.ERROR
    
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def parse_args(Args=None):
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Summarise multiple bigwig files into a single output bigwig file over a specified genomic region or whole genome.',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output bigwig file path'
    )
    
    parser.add_argument(
        '-i', '--input',
        nargs='+',
        required=True,
        help='Input bigwig file paths'
    )

    parser.add_argument(
        '--region',
        help='Genomic region in format chr:start-end (e.g., chr2:197,623,164-198,410,427). If not provided, will process whole genome.',
        default=None
    )

    parser.add_argument(
        '--function',
        choices=['sum', 'mean', 'median'],
        default='sum',
        help='Function to use for summarising bigwig values (default: sum)'
    )
    
    parser.add_argument(
        '--chunk-size',
        type=int,
        default=1000000,
        help='Window size in bp for reading bigwigs in chunks to manage memory usage (default: 1000000)'
    )
    
    parser.add_argument(
        '--precision',
        type=int,
        default=5,
        help='Number of decimal places to round output values to conserve file size (default: 5)'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='count',
        default=0,
        dest='verbosity',
        help='Increase verbosity (-v for WARNING, -vv for INFO, -vvv for DEBUG)'
    )
    
    if Args is None:
        return parser.parse_args()
    else:
        return parser.parse_args(Args)


def SummariseBigwigsInWindow(bigwig_out_fn, bigwigs_in_fn_list, window=None, function="sum", chunk_size=1000000, precision=5):
    """
    Given output bigwig filepath, and list of input bigwig filepaths, summarise the input bigwigs with a function. 
    For example, sum the values in all the input files. bigwig chrom headers should match.
    
    Parameters:
    -----------
    bigwig_out_fn : str
        Output bigwig file path
    bigwigs_in_fn_list : list
        List of input bigwig file paths
    window : str or None
        Genomic region in format chr:start-end. If None, processes whole genome.
    function : str
        Function to use for summarising ('sum', 'mean', 'median')
    chunk_size : int
        Size of chunks to read data in bp (for memory management)
    precision : int
        Number of decimal places to round output values
    """
    # Check that input bigwig headers chrom list matches across all files
    bigwigs = [pyBigWig.open(bw) for bw in bigwigs_in_fn_list]
    chroms = bigwigs[0].chroms()
    for i, bw in enumerate(bigwigs):
        if bw.chroms() != chroms:
            raise Exception(f'bigwig chromosome header list should match. {bigwigs_in_fn_list[0]} chroms do not match {bigwigs_in_fn_list[i]}.')
    
    # Initialize output file
    bw_out = pyBigWig.open(bigwig_out_fn, "w")
    bw_out.addHeader(list(chroms.items()))
    
    # Determine regions to process
    if window:
        # Single region specified
        chrom, coords = window.split(':')
        coords_float = [int(i.replace(',', '')) for i in coords.split('-')]
        start, end = coords_float
        base_regions = [(chrom, start, end)]
        logging.info(f"Processing region {chrom}:{start:,}-{end:,} (length: {end-start:,} bp)")
    else:
        # Whole genome - process all chromosomes
        base_regions = [(chrom, 0, size) for chrom, size in chroms.items()]
        total_bp = sum(size for _, _, size in base_regions)
        logging.info(f"Processing whole genome: {len(base_regions)} chromosomes, {total_bp:,} bp total")
    
    # Break base regions into chunks based on chunk_size
    regions_to_process = []
    for chrom, start, end in base_regions:
        region_length = end - start
        if region_length <= chunk_size:
            # Region is smaller than chunk size, process as-is
            regions_to_process.append((chrom, start, end))
            logging.debug(f"Added region {chrom}:{start:,}-{end:,} (length: {region_length:,} bp)")
        else:
            # Break region into chunks
            num_chunks = (region_length + chunk_size - 1) // chunk_size  # Ceiling division
            logging.info(f"Breaking {chrom}:{start:,}-{end:,} into {num_chunks} chunks of ~{chunk_size:,} bp")
            
            for chunk_i in range(num_chunks):
                chunk_start = start + (chunk_i * chunk_size)
                chunk_end = min(chunk_start + chunk_size, end)
                regions_to_process.append((chrom, chunk_start, chunk_end))
                logging.debug(f"Added chunk {chrom}:{chunk_start:,}-{chunk_end:,} (length: {chunk_end-chunk_start:,} bp)")
    
    logging.info(f"Total regions to process: {len(regions_to_process)}")
    logging.info(f"Using function: {function}")
    logging.info(f"Chunk size: {chunk_size:,} bp")
    logging.info(f"Output precision: {precision} decimal places")
    
    # Process each region (these are now already chunked appropriately)
    for region_i, (chrom, start, end) in enumerate(regions_to_process):
        logging.info(f"Processing region {region_i+1}/{len(regions_to_process)}: {chrom}:{start:,}-{end:,} (length: {end-start:,} bp)")
        
        # Collect arrays for this region
        region_arrays = []
        for bw in bigwigs:
            # Check if this bigwig has data in this region
            try:
                if bw.stats(chrom, start, end, type="max") == [None]:
                    logging.debug(f"No data in {bw} for region {chrom}:{start:,}-{end:,}")
                    continue
                else:
                    values = bw.values(chrom, start, end, numpy=True)
                    if values is not None:
                        region_arrays.append(values)
            except RuntimeError as e:
                # Handle cases where chromosome doesn't exist in bigwig
                logging.debug(f"Error accessing {chrom} in {bw}: {e}")
                continue
        
        if region_arrays:
            logging.debug(f"Writing region {chrom}:{start:,}-{end:,}")
            
            # Apply the specified function
            if function == "sum":
                region_summarised = numpy.sum(numpy.array(region_arrays), axis=0)
            elif function == "mean":
                region_summarised = numpy.mean(numpy.array(region_arrays), axis=0)
            elif function == "median":
                region_summarised = numpy.median(numpy.array(region_arrays), axis=0)
            else:
                raise ValueError(f"Unknown function: {function}")
            
            # Handle NaN values (replace with 0)
            region_summarised = numpy.nan_to_num(region_summarised, nan=0.0)
            
            # Round to specified precision to conserve file size
            region_summarised = numpy.round(region_summarised, decimals=precision)
            
            # Write to output bigwig
            bw_out.addEntries(chrom, start, values=region_summarised, span=1, step=1)
        else:
            logging.debug(f"No data found for region {chrom}:{start:,}-{end:,}")
    
    # Close files
    for bw in bigwigs:
        bw.close()
    bw_out.close()
    
    logging.info(f"Finished writing {bigwig_out_fn}")


def main(**kwargs):
    """Main function to execute bigwig summarisation."""
    SummariseBigwigsInWindow(
        bigwig_out_fn=kwargs['output'],  # Changed from 'output_bigwig'
        bigwigs_in_fn_list=kwargs['input'],  # Changed from 'input_bigwigs'
        window=kwargs['region'],
        function=kwargs['function'],
        chunk_size=kwargs['chunk_size'],
        precision=kwargs['precision']
    )


if __name__ == "__main__":
    # I like to script and debug with an interactive interpreter in the same
    # working dir as this script.. If using interactive interpreter, can step
    # thru the script with args defined below
    try:
        if sys.ps1:
            # Test basic sum functionality
            Args = "-vv --function mean -o scratch/test.bw -i test_data/Example_HTT_data/bigwigs/1_Fibroblast_polyA_A10_NA_1.bw test_data/Example_HTT_data/bigwigs/1_Fibroblast_polyA_A10_NA_2.bw test_data/Example_HTT_data/bigwigs/1_Fibroblast_polyA_A10_NA_3.bw".split(' ')
            # Test with region and mean function
            # Args = "-o scratch/test_mean.bw --region chr2:197,623,164-198,410,427 --function mean --chunk-size 500000 -vv -i bigwigs/H3K4ME3/NA18489.1.bw bigwigs/H3K4ME3/NA18498.2.bw".split(' ')
            # Test with median function and custom precision
            # Args = "-o scratch/test_median.bw --region chr2:197,623,164-198,410,427 --function median --precision 3 -vvv -i bigwigs/H3K4ME3/NA18489.1.bw bigwigs/H3K4ME3/NA18498.2.bw".split(' ')
            # Test whole genome with high precision
            # Args = "-o scratch/test_genome.bw --precision 8 --function sum -vv -i bigwigs/H3K4ME3/NA18489.1.bw bigwigs/H3K4ME3/NA18498.2.bw".split(' ')
            args = parse_args(Args=Args)
    except:
        args = parse_args()
    
    setup_logging(args.verbosity)
    main(**vars(args))