#! /usr/bin/env python3

import argparse

import misc.gffread_utils as gffread_utils
import misc.utils as utils

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script extracts all of the UTR information from a "
        "transcripts file created by gffread. It writes it into a csv file, so "
        "it can easily be imported and manipulated with pandas.")

    parser.add_argument('transcripts', help="A gffread-like fasta file")
    parser.add_argument('out', help="The (csv.gz) output file")


    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    utr_info = gffread_utils.get_all_utrs(args.transcripts)

    utils.write_df(utr_info, args.out, index=False)

if __name__ == '__main__':
    main()
