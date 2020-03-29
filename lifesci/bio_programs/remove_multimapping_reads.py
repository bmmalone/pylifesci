#! /usr/bin/env python3

import argparse
import lifesci.bam_utils as bam_utils

import logging
import pyllars.logging_utils as logging_utils
logger = logging.getLogger(__name__)

default_tmp = None

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script removes all reads with multiple alignments from a bam file. "
        "It then sorts and indexes the reads.")
    parser.add_argument('align_in', help="The input alignment file")
    parser.add_argument('align_out', help="The output alignment file with multimappers removed")
    parser.add_argument('--tmp', help= "The path where temporary files for samtools sort will "
        "be stored. If not given, then the samtools default tmp choice will be used.", 
        default=default_tmp)

    parser.add_argument('--do-not-call', action='store_true')

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    call = not args.do_not_call
    bam_utils.remove_multimapping_reads(args.align_in, args.align_out, call=call, tmp=args.tmp)

if __name__ == '__main__':
    main()
