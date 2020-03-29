#! /usr/bin/env python3

import argparse
import numpy as np

import pybedtools

import lifesci.bio as bio
import lifesci.bed_utils as bed_utils
import misc.parallel as parallel
import misc.utils as utils

default_num_cpus = 1

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script calculates the overlap, in base pairs, of two bed files. "
        "It also calculates the number of base pairs which are unique to each bed file.")

    parser.add_argument('bed_1', help="The first BED12+ file")
    parser.add_argument('bed_2', help="The second BED12+ file")

    parser.add_argument('--num-cpus', help="The number of CPUs to use when counting the "
        "BED12 feature lengths", type=int, default=default_num_cpus)
    
    args = parser.parse_args()

    programs = ['subtractBed']
    utils.check_programs_exist(programs)
        
    # read in the files and convert for use in bedtools
    bed_1_df = bed_utils.read_bed(args.bed_1)
    bed_2_df = bed_utils.read_bed(args.bed_2)

    bed_1 = pybedtools.BedTool.from_dataframe(bed_1_df[bed_utils.bed12_field_names])
    bed_2 = pybedtools.BedTool.from_dataframe(bed_2_df[bed_utils.bed12_field_names])

    # first, the bases unique to bed_1
    bed_1_only = bed_1.subtract(bed_2, split=True)
    bed_1_only_df = bed_1_only.to_dataframe(names=bed_utils.bed12_field_names)
    bed_1_only_sizes = parallel.apply_parallel(bed_1_only_df, args.num_cpus, bed_utils.get_bed_12_feature_length)

    bed_1_only_coverage = np.sum(bed_1_only_sizes)

    # now the bed_2 unique bases
    bed_2_only = bed_2.subtract(bed_1, split=True)
    bed_2_only_df = bed_2_only.to_dataframe(names=bio.bed12_field_names)
    bed_2_only_sizes = parallel.apply_parallel(bed_2_only_df, args.num_cpus, bed_utils.get_bed_12_feature_length)

    bed_2_only_coverage = np.sum(bed_2_only_sizes)

    # and the overlap
    bed_1_and_2 = bed_1.intersect(bed_2, split=True)
    bed_1_and_2_df = bed_1_and_2.to_dataframe(names=bed_utils.bed12_field_names)
    bed_1_and_2_sizes = parallel.apply_parallel(bed_1_and_2_df, args.num_cpus, bed_utils.get_bed_12_feature_length)

    bed_1_and_2_coverage = np.sum(bed_1_and_2_sizes)

    print("{} unique bases: {}".format(args.bed_1, bed_1_only_coverage))
    print("{} unique bases: {}".format(args.bed_2, bed_2_only_coverage))
    print("Overlapping bases: {}".format(bed_1_and_2_coverage))

if __name__ == '__main__':
    main()
