#! /usr/bin/env python3

import argparse
import logging
import pandas as pd

import lifesci.bio as bio
import lifesci.bed_utils as bed_utils
import misc.parallel as parallel
import pyllars.logging_utils as logging_utils
import misc.utils as utils

logger = logging.getLogger(__name__)

default_ccds_status_to_ignore = [
    'Withdrawn',
    'Reviewed, update pending'
]

default_num_cpus = 1

def parse_cds_locations(row, args):
    # make sure this is still a valid CCDS
    if row['ccds_status'] in args.ignore:
        return None
    
    # first, strip off the '[' and ']'
    cds_locations = row['cds_locations']
    cds_locations = cds_locations[1:-1]
    
    try:
        bed_exons = bed_utils.convert_genomic_coords_to_bed_blocks(cds_locations)
    except ValueError as ex:
        print(ex)
        print(row)
        bed_exons = None
        
    # and add the identifier so we can easily join later
    cds_id = row['gene'] + ":" + row['ccds_id']
    bed_exons['id'] = cds_id
    return bed_exons

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script converts the CCDS text files distributed by NCBI to "
        "valid BED12 files for use with other programs. Please see "
        "ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/README, \"CCDS.[YearMonthDay].txt\" for "
        "more information.")
    parser.add_argument('ccds', help="The CCDS.txt file downloaded from NCBI")
    parser.add_argument('out', help="The output bed.gz file")

    parser.add_argument('-i', '--ignore', help="The ccds_status entries to ignore.",
        default=default_ccds_status_to_ignore, nargs='*')

    parser.add_argument('-p', '--num-cpus', help="The number of processors to use "
        "for extracting the exon information", type=int, default=default_num_cpus)
    
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "Reading CCDS file"
    logger.info(msg)

    ccds_df = pd.read_csv(args.ccds, sep='\t')

    msg = "Copying simple values to BED"
    logger.info(msg)

    ccds_bed = pd.DataFrame()
    ccds_bed['seqname'] = ccds_df['#chromosome']
    ccds_bed['start'] = ccds_df['cds_from']
    ccds_bed['end'] = ccds_df['cds_to']
    ccds_bed['id'] = ccds_df['gene'] + ":" + ccds_df['ccds_id']
    ccds_bed['score'] = 0
    ccds_bed['strand'] = ccds_df['cds_strand']
    ccds_bed['thick_start'] = ccds_df['cds_from']
    ccds_bed['thick_end'] = ccds_df['cds_to']
    ccds_bed['color'] = 0

    msg = "Converting CCDS exons into BED12 blocks"
    logger.info(msg)

    cds_exon_info = parallel.apply_parallel(ccds_df, args.num_cpus, 
        parse_cds_locations, args, 
        progress_bar=True)

    cds_exon_info = [cei for cei in cds_exon_info if cei is not None]
    cds_exon_df = pd.DataFrame(cds_exon_info)

    msg = "Merging simple values and blocks"
    logger.info(msg)

    ccds_bed = ccds_bed.merge(cds_exon_df, on='id')

    # put the columns in the correct order
    ccds_bed = ccds_bed[bio.bed12_field_names]

    msg = "Writing the BED file"
    logger.info(msg)

    bed_utils.write_bed(ccds_bed, args.out)

if __name__ == '__main__':
    main()
