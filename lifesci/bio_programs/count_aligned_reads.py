#! /usr/bin/env python3

import argparse

import pandas as pd

import lifesci.bam_utils as bam_utils
import pyllars.pandas_utils as pd_utils
import pyllars.utils

import logging
import pyllars.logging_utils as logging_utils
logger = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script counts the number of uniquely- and multi-mapping "
        "reads in a list of bam files.")

    parser.add_argument('files', help="A glob-style re giving the filenames", 
        nargs='+')

    parser.add_argument('-o', '--out', help="A (csv.gz) output file", required=True)
    
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    datasets= []
    aligned_reads = []
    uniquely_aligned_reads = []

    for f in args.files:
        msg = "Processing file: {}".format(f)
        logger.info(msg)

        num_aligned_reads = bam_utils.count_aligned_reads(f)
        num_uniquely_aligned_reads = bam_utils.count_uniquely_mapping_reads(f)

        datasets.append(pyllars.utils.get_basename(f))
        aligned_reads.append(num_aligned_reads)
        uniquely_aligned_reads.append(num_uniquely_aligned_reads)

    msg = "Constructing data frame"
    logger.info(msg)

    df = pd.DataFrame()
    df['dataset'] = datasets
    df['aligned_reads'] = aligned_reads
    df['uniquely_aligned_reads'] = uniquely_aligned_reads
    df['multimapping_reads'] = df['aligned_reads'] - df['uniquely_aligned_reads']

    msg = "Writing data frame to disk"
    logger.info(msg)

    pd_utils.write_df(df, args.out, index=False)

if __name__ == '__main__':
    main()
