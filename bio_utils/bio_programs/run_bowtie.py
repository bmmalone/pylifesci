#! /usr/bin/env python3

import argparse
import os

import misc.bio as bio
import misc.shell_utils as shell_utils
import misc.slurm as slurm
import misc.utils as utils

import logging
import misc.logging_utils as logging_utils
logger = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script runs bowtie2 on all of the provided input files "
        "using the given index. By default, it does not save the alignments, "
        "aligned reads or unaligned reads. The respective flags must be given "
        "to retain the desired entities.")

    parser.add_argument('index', help="The bowtie2 index")
    parser.add_argument('out', help="The output directory")
    parser.add_argument('fastq', help="The fastq files", nargs='+')

    parser.add_argument('-a', '--alignments', help="If this flag is present, "
        "the alignments will be present in the output folder", action='store_true')
    parser.add_argument('--un-gz', help="If this flag is present, then the "
        "unaligned reads will be present in the output folder", action='store_true')
    parser.add_argument('--al-gz', help="If this flag is present, then the "
        "aligned reads will be present in the output folder", action='store_true')
    
    slurm.add_sbatch_options(parser)
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    programs = ['bowtie2', 'call-program']
    shell_utils.check_programs_exist(programs)

    if not os.path.exists(args.out):
        if not args.do_not_call:
            msg = "Creating output directory: {}".format(args.out)
            logger.info(msg)
            os.makedirs(args.out)

    for fastq in args.fastq:

        basename = utils.get_basename(fastq)

        out_files = []

        out = utils.abspath("dev","null") # we do not care about the alignments
        out_str = "-S {}".format(out)
        if args.alignments:
            n = "{}.bam".format(basename)
            out = os.path.join(args.out, n)
            out_str = "-S {}".format(out)
            out_files.append(out)

        un_gz_str = ""
        if args.un_gz:
            n = "{}.un-al.fastq.gz".format(basename)
            n = os.path.join(args.out, n)
            un_gz_str = "--un-gz {}".format(n)
            out_files.append(n)

        al_gz_str = ""
        if args.al_gz:
            n = "{}.al.fastq.gz".format(basename)
            n = os.path.join(args.out, n)
            al_gz_str = "--al-gz {}".format(n)
            out_files.append(n)

        cmd = "call-program bowtie2 -p {} --very-fast -x {} -U {} {} {} {}".format(
            args.num_cpus, args.index, fastq, out_str, un_gz_str, al_gz_str)

        slurm.check_sbatch(cmd, args=args)
        #shell_utils.check_call(cmd)

if __name__ == '__main__':
    main()
