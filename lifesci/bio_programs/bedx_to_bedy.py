#! /usr/bin/env python3

import argparse

import lifesci.bed_utils as bed_utils

import logging
import pyllars.logging_utils as logging_utils
logger = logging.getLogger(__name__)

default_bed_y = 12

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="The script \"downcasts\" a bed12+ file to a bed12 file. "
        "Optionally, it can down cast to fewer columns. By default, it also "
        "names the columns (using bed comments) according to lifesci."
        "bed_utils.bed12_field_names.")

    parser.add_argument('bed', help="The bed12+ file")
    parser.add_argument('out', help="The output bed12 (or less) file")

    parser.add_argument('--do-not-rename', help="If this flag is given, then "
        "the columns will not be renamed", action='store_true')

    parser.add_argument('-y', '--bed-y', help="The first y columns of the bed "
        "file will be retained, so this could also be used to create, for "
        "example, a bed6 from a bed12", type=int, default=default_bed_y)

    parser.add_argument('--do-not-compress', help="If this flag is given, then "
        "the output bed file will not be compressed. Either way this *does not* "
        "change the filename. If the file is to be compressed out should already "
        "include \".gz\".", action='store_true')
    
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    # get rid of the double negatives...
    rename = not args.do_not_rename
    compress = not args.do_not_compress

    msg = "Reading original bed"
    logger.info(msg)

    bed = bed_utils.read_bed(args.bed)

    msg = "Selecting columns to keep"
    logger.info(msg)

    num_columns_to_keep = args.bed_y
    if len(bed) < num_columns_to_keep:
        num_columns_to_keep = len(bed)

    if num_columns_to_keep > 12:
        msg = "bedx-to-bedy can only be used to keep at most 12 columns."
        raise ValueError(msg)

    bed = bed.iloc[:,:num_columns_to_keep]

    if rename:
        msg = "Renaming columns"
        logger.info(msg)
        bed.columns = bed_utils.bed12_field_names[:num_columns_to_keep]

    msg = "Writing \"down-cast\" bed to disk"
    logger.info(msg)

    bed_utils.write_bed(bed, args.out, compress=compress)

if __name__ == '__main__':
    main()
