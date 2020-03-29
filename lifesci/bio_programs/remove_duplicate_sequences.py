#! /usr/bin/env python3

import argparse
import logging
import lifesci.fastx_utils as fastx_utils

import pyllars.logging_utils as logging_utils

logger = logging.getLogger(__name__)

default_lower_precedence_re = None

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script removes all duplicate sequences from a list of fasta "
        "files and writes the remaining sequences back out as a fasta file."
        "\n\n"
        "If desired, a regular expression can be given for \"lower precedence\" "
        "sequence identifiers. An example of using this precedence operator "
        "is in removing duplicate sequences from a fasta file which combines "
        "de novo assembled transcripts and annotated ones. In case a de novo "
        "transcript matches an annotated one, we would prefer to keep only "
        "the annotated transcript and identifier. Thus, we would pass an RE "
        "matching the de novo assembled identifiers (which have a lower precedence)."
        "\n\n"
        "If a precedence re is not given, or two identifiers have the same "
        "precedence, the first identifier encountered will be kept.")

    parser.add_argument('fasta', help="The input fasta file(s)", nargs='+')
    parser.add_argument('-o', '--out', help="The output (fasta[.gz]) file",
        required=True)
    parser.add_argument('--compress', help="If this flag is given, then the output "
        "will be gzipped.", action='store_true')
    parser.add_argument('-l', '--lower-precedence-re', help="A regular expression "
        "that matches the identifiers of lower precendence transcripts. (See the "
        "description for more details.)", default=default_lower_precedence_re)

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    files = '\n'.join(args.fasta)
    msg = "Found the following files from the command line: {}".format(files)
    logger.info(msg)

    fastx_utils.remove_duplicate_sequences(args.fasta, args.out, compress=args.compress, 
        lower_precedence_re=args.lower_precedence_re, progress_bar=True)


if __name__ == '__main__':
    main()
