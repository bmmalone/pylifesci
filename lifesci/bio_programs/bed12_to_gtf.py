#! /usr/bin/env python3

import argparse

import logging
import misc.logging_utils as logging_utils
logger = logging.getLogger(__name__)

import pandas as pd

import lifesci.bed_utils as bed_utils
import lifesci.gtf_utils as gtf_utils

import misc.parallel as parallel

default_source = None
default_num_cpus = 1

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Convert a bed12 file into an equivalent gtf file. In "
        "particular, it uses the \"thick_start\" and \"thick_end\" fields to "
        "determine the CDS gtf entries. It only creates \"exon\" and \"CDS\" "
        "gtf entries. The bed columns after the standard 12 are included as "
        "attributes in the gtf file.")

    parser.add_argument('bed', help="The bed12 file. It must conform to the "
        "style expected by lifesci.bed_utils.")
    parser.add_argument('out', help="The (output) gtf file. It will conform "
        "to the style dictated by lifesci.gtf_utils.")

    parser.add_argument('-s', '--source', help="The name to use for the "
        "\"source\" column in the gtf file", default=default_source)

    parser.add_argument('-p', '--num-cpus', help="The number of CPUs to use "
        "for conversion", type=int, default=default_num_cpus)

    parser.add_argument('--add-gene-id', help="If this flag is present, then "
        "the \"id\" field will also be used as the \"gene_id\"",
        action='store_true')
    
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "Reading bed file"
    logger.info(msg)
    bed = bed_utils.read_bed(args.bed)

    if args.add_gene_id:
        msg = "Adding gene_id"
        logger.info(msg)
        bed['gene_id'] = bed['id']

    msg = "Expanding bed entries to gtf entries"
    logger.info(msg)

    gtf_entries = parallel.apply_parallel(
        bed,
        args.num_cpus,
        gtf_utils.get_gtf_entries,
        args.source,
        progress_bar = True
    )

    msg = "Joining gtf entries into large data frame"
    logger.info(msg)

    gtf_entries = pd.concat(gtf_entries)

    msg = "Sorting gtf entries"
    logger.info(msg)

    gtf_entries = gtf_entries.sort_values(['seqname', 'start', 'end'])
    gtf_entries = gtf_entries.reset_index(drop=True)

    msg = "Writing gtf to disk"
    logger.info(msg)
    gtf_utils.write_gtf(gtf_entries, args.out, compress=False)


if __name__ == '__main__':
    main()
