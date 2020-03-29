#! /usr/bin/env python3

import argparse
import logging

import pandas as pd

import lifesci.bam_utils as bam_utils
import lifesci.fastx_utils as fastx_utils
import misc.parallel as parallel

import pyllars.pandas_utils as pd_utils
import pyllars.logging_utils as logging_utils
import pyllars.utils

logger = logging.getLogger(__name__)

default_num_cpus = 1

file_type_choices = [
    "AUTO",
    "bam",
    "fasta",
    "fastq"
]

default_file_type = "AUTO"

auto_extensions = {
    "bam": set(["bam", "sam"]),
    "fastq": set(["fastq", "fastq.gz", "fq", "fq.gz"])
}

def guess_file_type(filename):
    file_type = "fasta"
    for ft, extensions in auto_extensions.items():
        for extension in extensions:
            if filename.endswith(extension):
                file_type = ft
                break
    return file_type

# file type-specific calls to the relevant length distribution calculators
def get_bam_length_distribution_df(fn):
    length_distribution_df = bam_utils.get_length_distribution(fn)
    basename = pyllars.utils.get_basename(fn)
    length_distribution_df['basename'] = basename
    return length_distribution_df

def get_fasta_length_distribution_df(fn):
    length_distribution_df = fastx_utils.get_length_distribution(fn, is_fasta=True)
    basename = pyllars.utils.get_basename(fn)
    length_distribution_df['basename'] = basename
    return length_distribution_df


def get_fastq_length_distribution_df(fn):
    length_distribution_df = fastx_utils.get_length_distribution(fn, is_fasta=False)
    basename = pyllars.utils.get_basename(fn)
    length_distribution_df['basename'] = basename
    return length_distribution_df

# function pointers to the appropriate length distribution calculator for each
# file type
file_type_get_length_distribution = {
    "bam": get_bam_length_distribution_df,
    "fasta": get_fasta_length_distribution_df,
    "fastq": get_fastq_length_distribution_df
}

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script counts the number of unique reads in the "
        "given files, all of which must be the same type. In the case of bam "
        "files, it only counts primary alignments (so it does not "
        "double-count multimappers, and it does not include unmapped reads "
        "present in the file.")

    parser.add_argument('files', help="The fasta, fastq or bam files",
        nargs='+')
    parser.add_argument('-o', '--out', help="The (csv.gz) output file "
        "containing the lengths and counts", required=True)

    parser.add_argument('-f', '--file-type', help="The type of the files. All "
        "files must be of the same type. If the \"AUTO\" file type is given, "
        "then the type will be guessed on the extension of the first file "
        "using the following heuristic: \"bam\" if the extension is\".bam\" "
        "or \".sam\"; "
        "\"fastq\" if the extension is \"fastq\", \"fastq.gz\", \"fq\", or "
        "\"fq.gz\"; \"fasta\" otherwise.", choices=file_type_choices, 
        default=default_file_type)

    parser.add_argument('-p', '--num-cpus', help="The number of CPUs to use",
        type=int, default=default_num_cpus)

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    if args.file_type == "AUTO":
        args.file_type = guess_file_type(args.files[0])
        msg = "The guessed file type is: {}".format(args.file_type)
        logger.info(msg)

    # grab the correct function pointer
    get_length_distribution = file_type_get_length_distribution[args.file_type]

    msg = "Collecting all read length distributions"
    logger.info(msg)

    all_length_distribution_dfs = parallel.apply_parallel_iter(
        args.files,
        args.num_cpus,
        get_length_distribution,
        progress_bar=True
    )

    msg = "Combining data frames into one large df"
    logger.info(msg)
    length_distribution_df = pd.concat(all_length_distribution_dfs)

    msg = "Writing counts to disk"
    logger.info(msg)

    pd_utils.write_df(length_distribution_df, args.out, index=False)

if __name__ == '__main__':
    main()
