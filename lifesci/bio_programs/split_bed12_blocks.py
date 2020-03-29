#! /usr/bin/env python3

import argparse
import logging
import pandas as pd

import lifesci.bio as bio
import lifesci.bed_utils as bed_utils

import pyllars.collection_utils as collection_utils
import pyllars.logging_utils as logging_utils

import misc.parallel as parallel

logger = logging.getLogger(__name__)

default_num_cpus = 1
default_num_groups = 500

def split_all_blocks(bed):
    exons = parallel.apply_df_simple(bed, bed_utils.split_bed12_blocks)
    exons = collection_utils.flatten_lists(exons)
    exons = pd.DataFrame(exons)
    return exons


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script splits BED12+ files into a BED6+ file. Each block "
        "(i.e., exon) in the original file is an individual feature in the new file. "
        "There are two extra fields, exon_index and transcript_start, which give the "
        "index of the exon within its transcript and the start of the exon in the "
        "\"spliced\" version of the transcript. The \"id\" column in the original file "
        "is used as the \"id\" in the new file, so the exons can easily be grouped.")

    parser.add_argument('bed', help="The BED12+ file")
    parser.add_argument('out', help="The output BED6+2 file")

    parser.add_argument('-p', '--num-cpus', help="The number of CPUs to use",
        type=int, default=default_num_cpus)

    parser.add_argument('--num-groups', help="The number of groups to split the "
        "bed file into for parallelization", type=int, default=default_num_groups)
    
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)
    
    msg = "Reading BED12+ file"
    logger.info(msg)

    bed = bed_utils.read_bed(args.bed)

    msg = "Splitting blocks"
    logger.info(msg)

    exons = parallel.apply_parallel_split(
        bed,
        args.num_cpus,
        #bio.split_bed12_blocks,
        split_all_blocks,
        progress_bar=True,
        num_groups = args.num_groups
    )

    msg = "Merging exons into a data frame"
    logger.info(msg)

    #exons = utils.flatten_lists(exons)
    #exons = pd.DataFrame(exons)
    exons = pd.concat(exons)

    fields = bed_utils.bed6_field_names + ['exon_index', 'transcript_start']
    exons = exons[fields]

    msg = "Writing BED6+2 file"
    logger.info(msg)

    bed_utils.write_bed(exons, args.out)

if __name__ == '__main__':
    main()
