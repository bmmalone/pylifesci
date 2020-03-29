#! /usr/bin/env python3

import argparse
from contextlib import contextmanager, closing
import os
import sys
import tempfile

import pysam

import misc.parallel as parallel
import misc.shell_utils as shell_utils
import misc.slurm as slurm
import misc.utils as utils


import logging
import misc.logging_utils as logging_utils
logger = logging.getLogger(__name__)

default_chrom = 'all'
default_start = 0
default_end = None


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script converts bam files to bigWig files. It is mostly "
        "a port of this script: https://github.com/chapmanb/bcbb/blob/master/nextgen/scripts/bam_to_wiggle.py "
        "by Brad Chapman which avoids a few dependencies.\n\nThe wigToBigWig "
        "program (from UCSC tools) must be in the path.\n\nN.B. If given, the "
        "start and end coordinates must be base-0.")

    parser.add_argument('bam', help="The bam file", nargs='+')
    parser.add_argument('-o', '--overwrite', help="If this flag is given, then "
        "the bigWig file will be created whether it exists or not", 
        action='store_true')
    parser.add_argument('-c', '--chrom', help="If specified, only alignments "
        "from this chromosome will be in the output", default=default_chrom)
    parser.add_argument('-s', '--start', help="If specied, only alignments "
        "from this position will be in the output", default=default_start)
    parser.add_argument('-e', '--end', help="If specied, only alignments "
        "up to this position will be in the output", default=default_end)

    parser.add_argument('-n', '--normalize', help="If this flag is given, "
        "then values will be normalized to reads per million", action='store_true')

    parser.add_argument('-t', '--use-tempfile', help="If this flag is given, "
        "then a temp file will be used to avoid permission issues", 
        action='store_true')

    parser.add_argument('-k', '--keep-wig', help="If this flag is given, then "
        "the wiggle file will not be deleted", action='store_true')
    
    slurm.add_sbatch_options(parser)
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    programs = ['wigToBigWig']
    shell_utils.check_programs_exist(programs)

    if args.use_slurm:
        cmd = ' '.join(sys.argv)
        slurm.check_sbatch(cmd, args=args)
        return

    parallel.apply_parallel_iter(args.bam, args.num_cpus, bam_to_wiggle, args,
        progress_bar=True)

def bam_to_wiggle(bam, args):
    out = "{}.bigWig".format(bam)
   
    regions = [(args.chrom, args.start, args.end)]

    bigWig_exists =  (os.path.exists(out) and os.path.getsize(out) > 0)
    if args.overwrite or not bigWig_exists:
        if args.use_tempfile:
            #Use a temp file to avoid any possiblity of not having write permission
            out_handle = tempfile.NamedTemporaryFile(delete=False)
            wig_file = out_handle.name
        else:
            out_base = os.path.splitext(out)[0]
            wig_file = "{}.wig".format(out_base)
            out_handle = open(wig_file, "w")

        msg = "Writing bam to wig"
        logger.info(msg)

        with closing(out_handle):
            chr_sizes, wig_valid = write_bam_track(bam, regions, out_handle,
                                                   args.normalize)
        try:
            msg = "Converting wig to bigWig"
            logger.info(msg)

            if wig_valid:
                convert_to_bigwig(wig_file, out, chr_sizes, args)
        finally:
            if not args.keep_wig:
                msg = "Removing wig file"
                logger.info(msg)

                os.remove(wig_file)

            else:
                msg = "Keeping wig file"
                logger.info(msg)

    else:
        msg = "The bigWig file already exists. Quitting."
        logger.warning(msg)

    
@contextmanager
def indexed_bam(bam_file):
    if not os.path.exists(bam_file + ".bai"):
        pysam.index(bam_file)
    sam_reader = pysam.Samfile(bam_file, "rb")
    yield sam_reader
    sam_reader.close()


def write_bam_track(bam_file, regions, out_handle, normalize):
    out_handle.write("track %s\n" % " ".join(["type=wiggle_0",
        "name=%s" % os.path.splitext(os.path.split(bam_file)[-1])[0],
        "visibility=full",
        ]))
    normal_scale = 1e6
    is_valid = False
    with indexed_bam(bam_file) as work_bam:
        total = sum(1 for r in work_bam.fetch() if not r.is_unmapped) if normalize else None
        sizes = list(zip(work_bam.references, work_bam.lengths))
        if len(regions) == 1 and regions[0][0] == "all":
            regions = [(name, 0, length) for name, length in sizes]
        for chrom, start, end in regions:
            if end is None and chrom in work_bam.references:
                end = work_bam.lengths[work_bam.references.index(chrom)]
            assert end is not None, "Could not find %s in header" % chrom
            out_handle.write("variableStep chrom=%s\n" % chrom)
            for col in work_bam.pileup(chrom, start, end):
                if normalize:
                    n = float(col.n) / total * normal_scale
                else:
                    n = col.n
                out_handle.write("%s %.1f\n" % (col.pos+1, n))
                is_valid = True
    return sizes, is_valid

def convert_to_bigwig(wig_file, out, chr_sizes, args):
    size_file = "%s-sizes.txt" % (os.path.splitext(wig_file)[0])
    with open(size_file, "w") as out_handle:
        for chrom, size in chr_sizes:
            out_handle.write("%s\t%s\n" % (chrom, size))
    try:
        msg = "Calling wigToBigWig"
        logger.info(msg)

        cmd = "wigToBigWig {} {} {}".format(wig_file, size_file, out)
        utils.check_call(cmd)
    finally:
        if args.keep_wig:
            msg = "Keeping size file"
            logger.info(msg)

        else:
            msg = "Removing size file"
            logger.info(msg)

            os.remove(size_file)
    return out

if __name__ == '__main__':
    main()
