#! /usr/bin/env python3

import argparse

import pandas as pd

import lifesci.fastx_utils as fastx_utils
import pyllars.utils
import pyllars.pandas_utils as pd_utils

import logging
import pyllars.logging_utils as logging_utils
logger = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script counts the number of reads in the given fastq "
        "(possibly gzipped) files.")

    parser.add_argument('files', help="A glob-style re giving the filenames", 
        nargs='+')

    parser.add_argument('-o', '--out', help="A (csv.gz) output file", required=True)
    
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    datasets = []
    read_counts = []

    for f in args.files:
        msg = "Processing file: {}".format(f)
        logger.info(msg)

        read_count = fastx_utils.get_read_count(f, is_fasta=False)

        datasets.append(pyllars.utils.get_basename(f))
        read_counts.append(read_count)

    df = pd.DataFrame()
    df['dataset'] = datasets
    df['reads'] = read_counts

    msg = "Writing data frame to disk"
    logger.info(msg)

    pd_utils.write_df(df, args.out, index=False)


if __name__ == '__main__':
    main()
