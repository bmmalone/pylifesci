#! /usr/bin/env python3

import argparse
import os
import pandas as pd

import lifesci.bed_utils as bed_utils
import pyllars.collection_utils as collection_utils

import misc.parallel as parallel

import logging
import pyllars.logging_utils as logging_utils
logger = logging.getLogger(__name__)

default_min_a_overlap = 0
default_min_b_overlap = 0
default_num_cpus = 1

def split_all_blocks(bed):
    exons = parallel.apply_df_simple(bed, bed_utils.split_bed12_blocks)
    exons = collection_utils.flatten_lists(exons)
    exons = pd.DataFrame(exons)
    return exons


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script removes all of the entries from A which overlap "
        "any of the entries from B. Optionally, some minimum overlap can be "
        "given, in terms of overlap fraction.")

    parser.add_argument('bed_a', help="The bed file from which entries will be "
        "removed")
    parser.add_argument('bed_b', help="The bed file used to find entries in A "
        "to remove")

    parser.add_argument('out', help="The output (bed.gz) file")

    parser.add_argument('--min-a-overlap', help="A minimum fraction required "
        "for overlap of the A entries. \"0\" means \"at least one bp.\"",
        type=float, default=default_min_a_overlap)
    
    parser.add_argument('--min-b-overlap', help="A minimum fraction required "
        "for overlap of the B entries. \"0\" means \"at least one bp.\"",
        type=float, default=default_min_b_overlap)

    parser.add_argument('--split', help="If this flag is given, then the bed "
        "entries in both files will be split. This can be somewhat slow, "
        "depending on the number of entries in the files.", action='store_true')

    parser.add_argument('--exons', help="If the bed entries have already been "
        "split and the exon bed6+2 file (from split-bed12-blocks program) is "
        "available, then that can be given with this option. The exons from "
        "that file will be used for both A and B.", default=None)

    parser.add_argument('--exons-a', help="As with the --exons argument, but "
        "these exons will only be used for A", default=None)

    parser.add_argument('--exons-b', help="As with the --exons argument, but "
        "these exons will only be used for B", default=None)

    parser.add_argument('-p', '--num-cpus', help="The number of CPUs to use for "
        "certain parts of the script", type=int, default=default_num_cpus)
    
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)
    
    exons_given = args.exons is not None
    exons_a_given = args.exons_a is not None
    exons_b_given = args.exons_b is not None

    # check if the exons files exist
    if exons_given and (not os.path.exists(args.exons)):
        msg = "The exons file does not exist: {}".format(args.exons)
        raise FileNotFoundError(msg)

    if exons_a_given and (not os.path.exists(args.exons_a)):
        msg = "The exons_a file does not exist: {}".format(args.exons_a)
        raise FileNotFoundError(msg)

    if exons_b_given and (not os.path.exists(args.exons_b)):
        msg = "The exons_b file does not exist: {}".format(args.exons_b)
        raise FileNotFoundError(msg)

    exons_a_only = exons_a_given and not exons_b_given
    exons_b_only = not exons_a_given and exons_b_given
    if exons_a_only or exons_b_only:
        msg = ("Only one of --exons-a, --exons-b was given. This is valid, "
            "but please ensure this is the desired behavior.")
        logger.warning(msg)

    # make sure we weren't given contradictory flags
    if args.split and exons_given:
        msg = "Both --split and --exons were given. Only one of these is allowed."
        raise ValueError(msg)

    if exons_given and (exons_a_given or exons_b_given):
        msg = ("Both --exons and (--exons-a or --exons-b) were given. --exons "
            "should not be given with the --exons-a and --exons-b arguments.")
        raise ValueError(msg)

    exons = None
    exons_a = None
    exons_b = None

    msg = "Reading bed A"
    logger.info(msg)
    bed_a = bed_utils.read_bed(args.bed_a)

    msg = "Reading bed B"
    logger.info(msg)
    bed_b = bed_utils.read_bed(args.bed_b)

    if args.split:
        msg = "Splitting bed A"
        logger.info(msg)

        exons_a = parallel.apply_parallel_split(
            bed_a,
            args.num_cpus,
            split_all_blocks,
            progress_bar=True,
            num_groups = args.num_groups
        )

        exons_a = pd.concat(exons_a)

        msg = "Splitting bed B"
        logger.info(msg)

        exons_b = parallel.apply_parallel_split(
            bed_b,
            args.num_cpus,
            split_all_blocks,
            progress_bar=True,
            num_groups = args.num_groups
        )

        exons_b = pd.concat(exons_b)

    if exons_given:
        msg = "Reading exons"
        logger.info(msg)
        exons = bed_utils.read_bed(args.exons)

    if exons_a_given:
        msg = "Reading A exons"
        logger.info(msg)
        exons_a = bed_utils.read_bed(args.exons_a)

    if exons_b_given:
        msg = "Reading B exons"
        logger.info(msg)
        exons_b = bed_utils.read_bed(args.exons_b)

    msg = "Finding all A entries which overlap B entries"
    logger.info(msg)

    remaining_a_ids = bed_utils.subtract_bed(
                    bed_a, 
                    bed_b, 
                    min_a_overlap=args.min_a_overlap, 
                    min_b_overlap=args.min_b_overlap, 
                    exons=exons,
                    exons_a=exons_a,
                    exons_b=exons_b)

    msg = "Filtering the A entries which had overlaps"
    logger.info(msg)

    m_remaining = bed_a['id'].isin(remaining_a_ids)
    bed_a_remaining = bed_a[m_remaining]

    msg = "Writing remaining A entries to disk"
    logger.info(msg)

    bed_utils.write_bed(bed_a_remaining, args.out)

if __name__ == '__main__':
    main()
