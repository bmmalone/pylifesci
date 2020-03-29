#! /usr/bin/env python3

import argparse
import logging

import lifesci.bed_utils as bed_utils
import lifesci.fastx_utils as fastx_utils
import pyllars.logging_utils as logging_utils

logger = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script extracts the (spliced) sequences from a bed file. "
        "Unless otherwise specified, it treats the input file as a bed12 file and "
        "extracts the spliced sequences. It reverse complements sequences from the "
        "reverse strand. The output is a fasta file. The headers in the fasta file "
        "are the \"id\"s (fourth column) from the bed file.")

    parser.add_argument('bed', help="The bed file")
    parser.add_argument('fasta', help="The chromosome fasta file")
    parser.add_argument('out', help="The output fasta file")

    parser.add_argument('--bed4', help="If this flag is present, the file will be "
        "treated as a bed4+ file. That is, the exon blocks will not be considered.",
        action='store_true')
    
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)
    
    # we will split the sequences unless the --bed4 flag was given
    split_exons = not args.bed4

    all_transcript_sequences = bed_utils.get_all_bed_sequences(args.bed, args.fasta, split_exons)

    msg = "Writing transcript sequences to disk"
    logger.info(msg)
    fastx_utils.write_fasta(all_transcript_sequences, args.out, compress=False)

if __name__ == '__main__':
    main()
