#! /usr/bin/env python3

import argparse
import pandas as pd

import lifesci.bed_utils as bed_utils

import logging
import pyllars.logging_utils as logging_utils
logger = logging.getLogger(__name__)


DUPLICATE_FIELDS = [
    'seqname', 
    'start', 
    'end', 
    'strand', 
    'num_exons', 
    'exon_lengths', 
    'exon_genomic_relative_starts'
]


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Remove duplicate entries from a list of bed12+ files. "
        "Write the non-redundant entries to a new file. Precedence among "
        "duplicates is arbitrary.")

    parser.add_argument('bed', help="The input bed file(s)", nargs='+')
    parser.add_argument('-o', '--out', help="The output bed(.gz) file",
        required=True)

    parser.add_argument('--compress', help="If this flag is given, the output "
        "will be gzipped. The output filename *will not* be changed (so it "
        "should already end in \".gz\").", action='store_true')
    
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "Reading bed files"
    logger.info(msg)

    all_bed = [
        bed_utils.read_bed(b) for b in args.bed
    ]

    for f, b in zip(args.bed, all_bed):
        msg = "{}. number of entities: {}".format(f, len(b))
        logger.debug(msg)

    msg = "Concatenating bed entries"
    logger.info(msg)
    all_bed_df = pd.concat(all_bed)

    msg = "Removing duplicate entries"
    logger.info(msg)
    all_bed_df = all_bed_df.drop_duplicates(subset=DUPLICATE_FIELDS)

    msg = "number of non-redundant entries: {}".format(len(all_bed_df))
    logger.debug(msg)

    msg = "Sorting non-redundant entries"
    logger.info(msg)
    sort_fields = ['seqname', 'start', 'end', 'strand']
    all_bed_df = all_bed_df.sort_values(by=sort_fields)

    msg = "Writing sorted, non-redundant entries to disk"
    logger.info(msg)
    bed_utils.write_bed(all_bed_df, args.out, compress=args.compress)

if __name__ == '__main__':
    main()
