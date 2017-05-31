#! /usr/bin/env python3

import argparse
import logging
import pandas as pd

import misc.bio as bio
import misc.parallel as parallel
import misc.utils as utils

default_max_size = 500e6
default_num_procs = 1

def get_original_coordinates(row, max_size):
    seqname = row['seqname']
    split = seqname.split("_")

    original_seqname = split[0]
    split_num = int(split[1])

    split_offset = split_num * max_size

    original_start = row['start'] + split_offset
    original_end = row['end'] + split_offset

    original_thick_start = row['thick_start'] + split_offset
    original_thick_end = row['thick_end'] + split_offset

    ret = {
        'seqname': original_seqname,
        'start': original_start,
        'end': original_end,
        'thick_start': original_thick_start,
        'thick_end': original_thick_end
    }

    return ret

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script rejoins annotations which were split with the "
        "split-long-chromosomes script. Importantly, the \"max-size\" parameter here "
        "must match the \"max-parameter\" used for that script.")
    
    parser.add_argument('bed', help="A BED file using the annotations from a GTF file "
        "created by split-long-chromosomes")
    parser.add_argument('out', help="A BED file in which the chromosomes have been rejoined")
    parser.add_argument('--max-size', help="The largest allowed size (in bp) for a "
        "chromosome", type=int, default=default_max_size)

    parser.add_argument('--num-procs', help="The number of processors to use", type=int,
        default=default_num_procs)
    
    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)
   
    msg = "Reading BED"
    logging.info(msg)

    bed = bio.read_bed(args.bed)

    msg = "Updating BED coordinates"
    logging.info(msg)

    original_coordinates = parallel.apply_parallel(bed, args.num_procs,
                                get_original_coordinates, args.max_size,
                                progress_bar=True)

    original_coordinates_df = pd.DataFrame(original_coordinates)

    bed['seqname'] = original_coordinates_df['seqname']
    bed['start'] = original_coordinates_df['start'].astype(int)
    bed['end'] = original_coordinates_df['end'].astype(int)
    bed['thick_start'] = original_coordinates_df['thick_start'].astype(int)
    bed['thick_end'] = original_coordinates_df['thick_end'].astype(int)

    msg = "Writing updated BED to disk"
    logging.info(msg)

    bio.write_bed(bed, args.out)

if __name__ == '__main__':
    main()
