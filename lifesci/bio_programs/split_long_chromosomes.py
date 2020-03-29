#! /usr/bin/env python3

from io import StringIO
from functools import partial
import logging

import numpy as np
import pandas as pd

import lifesci.fastx_utils as fastx_utils
import lifesci.gtf_utils as gtf_utils

import pyllars.logging_utils as logging_utils

import argparse

default_max_size = 500e6

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script splits large chromosomes into several smaller ones. This "
        "is required to use BedTools for chromosomes with size larger than about 500M. In "
        "particular, this script splits chromosomal sequences into smaller chunks and "
        "updates GTF annotations to use the smaller chromosomes.\n\nFor more information, "
        "see https://groups.google.com/forum/#!topic/bedtools-discuss/t-nQSCxaFGE")
    parser.add_argument('fasta', help="The chromosome sequence file")
    parser.add_argument('gtf', help="The annotation file")
    parser.add_argument('out', help="The base output files. The script will create the "
        "files <out>.fa and <out>.gtf.")
    parser.add_argument('--max-size', help="The largest allowed size (in bp) for a "
        "chromosome", type=int, default=default_max_size)
    
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "Splitting fasta sequences"
    logging.info(msg)

    fasta = fastx_utils.get_read_iterator(args.fasta, is_fasta=True)
    split_fasta = {}

    for name, seq in fasta:
        name_split = name.split(" ")
        split_seqs = [l for l in iter(partial(StringIO(seq).read, int(args.max_size)), '')]
        
        for i, split_seq in enumerate(split_seqs):
            n = name_split[0] # a bit of a hack to save the base sequence name
            name_split[0] = "{}_{}".format(name_split[0], i)
            split_name = ' '.join(name_split)
            
            split_fasta[split_name] = split_seq
            
            name_split[0] = n

    msg = "Writing fasta output file"
    logging.info(msg)

    fasta_out = "{}.fa".format(args.out)
    fastx_utils.write_fasta(split_fasta, fasta_out, compress=False, progress_bar=True)

    msg = "Reading GTF"
    logging.info(msg)

    gtf = gtf_utils.read_gtf(args.gtf)

    msg = "Updating GTF coordinates"
    logging.info(msg)

    # get the split for each feature
    start_split_num = gtf['start'] // args.max_size
    end_split_num = gtf['end'] // args.max_size
    
    # wrap the coordinates based on the max size
    gtf['start'] = np.mod(gtf['start'], args.max_size)
    gtf['end'] = np.mod(gtf['end'], args.max_size)

    # for the names, we need them as strings without ".0" at the end
    start_split_num = start_split_num.astype(int)
    start_split_num = start_split_num.astype(str)

    end_split_num = end_split_num.astype(int)
    end_split_num = end_split_num.astype(str)

    # create the start and end names
    split_start_seqname = gtf['seqname'] + "_" + start_split_num
    split_end_seqname = gtf['seqname'] + "_" + end_split_num

    gtf['start_seqname'] = split_start_seqname
    gtf['end_seqname'] = split_end_seqname

    # remove features which span the gaps
    m_span = gtf['start_seqname'] != gtf['end_seqname']
    gtf = gtf[~m_span]

    num_spanning = sum(m_span)
    msg = ("Number of features spanning length boundaries: {}.\n\nThese features "
        "will be discarded.".format(num_spanning))
    logging.warning(msg)

    # and update the seqnames
    gtf['seqname'] = gtf['start_seqname']

    msg = "Writing GTF output file"
    logging.info(msg)

    gtf_out = "{}.gtf".format(args.out)
    gtf_utils.write_gtf(gtf, gtf_out, compress=False)


if __name__ == '__main__':
    main()
