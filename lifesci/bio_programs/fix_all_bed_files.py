#! /usr/bin/env python3

import argparse
import glob

import lifesci.bed_utils as bed_utils

import logging
import pyllars.logging_utils as logging_utils
logger = logging.getLogger(__name__)

default_bed_extensions = ["bed", "bed.gz"]

def fix_bed(bed_file):
    msg = "(please type y/n) "
    print(msg)

    res = ""
    while res not in ["y", "n"]:
        res = input()

    if res == 'y':
        return True
    else:
        return False

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script replaces all score and color values in bed files with "
        "'0'. It applies this to all bed files in the current directory.")

    parser.add_argument('--no-ask', help="By default, the program will ask to replace "
        "the values for each bed file. If this flag is given, then the asking will be "
        "skipped.", action='store_true')

    parser.add_argument('--bed-extensions', help="The extensions to treat as "
        "bed files", nargs='+', default=default_bed_extensions)
    
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    ask = not args.no_ask

    for bed_extension in args.bed_extensions:
        re = "*{}".format(bed_extension)
        bed_files = glob.glob(re)

        for bed_file in bed_files:
            print("fix: {}".format(bed_file))
            if (not ask) or fix_bed(bed_file):
                bed = bed_utils.read_bed(bed_file)
                bed['score'] = 0
                bed['color'] = 0

                bed_utils.write_bed(bed, bed_file)

if __name__ == '__main__':
    main()
