#! /usr/bin/env python3

import argparse
from contextlib import ExitStack

import lifesci.bam_utils as bam_utils

import logging
import pyllars.logging_utils as logging_utils
logger = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script filters bam A by the read ids in bam B. In "
        "particular, only the reads in A with ids *which appear* in B are kept.")

    parser.add_argument('bam_a', help="The bam file to filter")
    parser.add_argument('bam_b', help="The bam file whose ids will be kept in A")
    parser.add_argument('bam_out', help="The output (bam) file")
    
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "Reading identifiers from B"
    logger.info(msg)
    ids_to_keep = bam_utils.get_read_identifiers(args.bam_b)

    msg = "Filtering reads from A which do not appear in B"
    logger.info(msg)

    with ExitStack() as stack:
        bam_a = stack.enter_context(bam_utils.get_pysam_alignment_file(args.bam_a))
        bam_out = stack.enter_context(bam_utils.get_pysam_alignment_file(
            args.bam_out, "wb", template=bam_a))

        for read in bam_a.fetch():
            if read.query_name in ids_to_keep:
                bam_out.write(read)

if __name__ == '__main__':
    main()
