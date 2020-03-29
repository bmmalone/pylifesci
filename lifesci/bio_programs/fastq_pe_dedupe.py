#! /usr/bin/env python3

import argparse

import collections
import lifesci.fastx_utils as fastx_utils
import pyllars.utils
import tqdm

from contextlib import ExitStack

import logging
import pyllars.logging_utils as logging_utils
logger = logging.getLogger(__name__)

stored_read_pair = collections.namedtuple(
    "stored_read_pair",
    "r1_name,r2_name,r1_qual,r2_qual"
)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script removes duplicate reads from paired-end fastq files. "
        "It only treats exact matches as duplicates (i.e., if reads are of different "
        "lengths, it does not consider exact substring matches as duplicates). The "
        "duplicate with the highest average quality score is retained.\n\nThis script "
        "is not designed to work with fasta files.")

    parser.add_argument('fastq_1', help="The first mate file")
    parser.add_argument('fastq_2', help="The second mate file")

    parser.add_argument('out_1', help="The de-duped first mate file")
    parser.add_argument('out_2', help="The de-duped second mate file")

    parser.add_argument('--do-not-compress', help="Unless this flag is given, the "
        "output will be gzipped. N.B. \".gz\" *will not* be adde to the file names.",
        action='store_true')
    
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    compress = not args.do_not_compress

    msg = "Counting reads in each file"
    logger.info(msg)

    # we will use the counts for displaying a progress bar
    num_reads_1 = fastx_utils.get_read_count(args.fastq_1, is_fasta=False)
    num_reads_2 = fastx_utils.get_read_count(args.fastq_2, is_fasta=False)

    # but avoid errors by making sure the counts match
    if num_reads_1 != num_reads_2:
        msg = "The number of reads in the files do not match ({} vs. {})".format(
            num_reads_1, num_reads_2)
        raise ValueError(msg)

    msg = "Creating read iterators"
    logger.info(msg)

    fastq_1_iter = fastx_utils.get_read_iterator(args.fastq_1, is_fasta=False)
    fastq_2_iter = fastx_utils.get_read_iterator(args.fastq_2, is_fasta=False)
    fastq_iter = zip(fastq_1_iter, fastq_2_iter)

    msg = "Detecting duplicates"
    logger.info(msg)

    seen_reads = {}

    for r1, r2 in tqdm.tqdm(fastq_iter, total=num_reads_1):
        r_key = (r1[1], r2[1])
        
        prev_val = seen_reads.get(r_key, None)
        
        if prev_val is None:
            srp = stored_read_pair(
                r1_name = r1[0],
                r2_name = r2[0],
                r1_qual = r1[2],
                r2_qual = r2[2]
            )
            seen_reads[r_key] = srp
        else:
            new_qual_score = sum(r1[2].encode()) + sum(r1[2].encode())
            prev_qual_score = (sum(prev_val.r1_qual.encode()) + 
                                    sum(prev_val.r2_qual.encode()))
            
            if new_qual_score > prev_qual_score:
                srp = stored_read_pair(
                    r1_name = r1[0],
                    r2_name = r2[0],
                    r1_qual = r1[2],
                    r2_qual = r2[2]
                )
                seen_reads[r_key] = srp

    msg = "Writing the de-duped files to disk"
    logger.info(msg)

    with ExitStack() as stack:
        out_1 = pyllars.utils.open_file(args.out_1, 'w', compress=compress)
        out_2 = pyllars.utils.open_file(args.out_2, 'w', compress=compress)
        
        for (seqs, srp) in tqdm.tqdm(seen_reads.items()):
            fastx_utils._write_fastq_entry(out_1, srp.r1_name, seqs[0], srp.r1_qual)
            fastx_utils._write_fastq_entry(out_2, srp.r2_name, seqs[1], srp.r2_qual)

if __name__ == '__main__':
    main()
