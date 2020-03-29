#! /usr/bin/env python3

import argparse
import logging
import textwrap

import tqdm

import misc.bio as bio
import misc.utils as utils

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script reorders sequences in a fasta file to be in the same "
        "order as the sequences in the STAR transcriptInfo.tab file.")
    parser.add_argument('transcript_info', help="The STAR transcriptInfo.txt file, which "
        "contains the desired sequence order")
    parser.add_argument('fasta', help="The fasta file containing the sequences. The index "
        "should already be created, and the keys for the index must exactly match the "
        "sequence names in the transcript_info file.")
    parser.add_argument('out', help="The output fasta file with the sequences rearranged")
    parser.add_argument('--compress', help="If this flag is present, then the output "
        "will be gzipped", action='store_true')

    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    msg = "Reading STAR transcript file"
    logging.info(msg)
    transcript_info = bio.read_star_tr_file(args.transcript_info)

    msg = "Reading transcript fasta file"
    logging.info(msg)
    fasta = bio.get_fasta_dict(args.fasta) #
    trans_id_mapping = {
        x.split()[0] : x
            for x in fasta.keys()
    }

    # open the tab file    
    if args.compress:
        out = gzip.open(args.out, 'wt')
    else:
        out = open(args.out, 'w')

    for i in tqdm.trange(len(transcript_info)):
        trans_id = transcript_info.iloc[i]['ID']
        fasta_id = trans_id_mapping[trans_id]

        header = ">{}\n".format(fasta_id)
        out.write(header)
        seq = str(fasta[fasta_id])

        wrapped_seq = textwrap.fill(seq)
        out.write(wrapped_seq)
        out.write("\n")

    out.close()

if __name__ == '__main__':
    main()
