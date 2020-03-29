#! /usr/bin/env python3

import argparse
from contextlib import ExitStack
import numpy as np

from Bio import SeqIO

import misc.utils as utils

import logging
import misc.logging_utils as logging_utils
logger = logging.getLogger(__name__)

default_max_quality_score = 40
default_quality_mean = 40
default_quality_std = 1

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script converts a fasta file to a fastq file with dummy "
        "quality scores. It is adapted from https://www.biostars.org/p/99886/")
    
    parser.add_argument('fasta', help="The input fasta file")
    parser.add_argument('out', help="The (output) fastq file")
    parser.add_argument('-c', '--compress', help="If this flag is given, then "
        "the output file will be compressed", action='store_true')

    parser.add_argument('-q', '--quality-mean', help="The mean quality score to "
        "use in the fastq file", type=int, default=default_quality_mean)

    parser.add_argument('--quality-std', help="The standard deviation in the "
        "quality scores", type=int, default=default_quality_std)

    parser.add_argument('--max-quality-score', help="The maximum quality score. "
        "Scores sampled higher than this will be floored.", type=int, 
        default=default_max_quality_score)

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    with ExitStack() as stack:
        fasta = stack.enter_context(open(args.fasta, "r"))
        fastq = stack.enter_context(utils.open(args.out, "w", compress=args.compress))

        for record in SeqIO.parse(fasta, 'fasta'):
            quals = np.random.normal(args.quality_mean, args.quality_std, len(record))
            quals = np.clip(quals, 0, 40)
            record.letter_annotations["phred_quality"] = quals
            SeqIO.write(record, fastq, "fastq")

if __name__ == '__main__':
    main()
