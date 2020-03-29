#! /usr/bin/env python3

import argparse
import misc.gffread_utils as gffread_utils
import misc.bio as bio
import pandas as pd     

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script extracts the CDS coordinates from the fasta sequence "
        "headers. It expects the format produces by gffread -W -w (from cufflinks).")
    parser.add_argument('fasta', help="The fasta file")
    parser.add_argument('out', help="The output file, in BED3 format")
    args = parser.parse_args()

    headers = gffread_utils.get_all_headers(args.fasta)

    # filter everything which does not have a cds_start and cds_end
    m_no_start = headers['cds_start'] == 0
    m_no_end= headers['cds_end'] == 0
    headers = headers[~(m_no_start & m_no_end)]

    # subtract 1 from start because the gff indices start at 1
    # no need to subtract from the end because bed is open at the end
    headers['cds_start'] = headers['cds_start'] - 1

    cds_only_df = pd.DataFrame()
    cds_only_df['seqname'] = headers['seqname']
    cds_only_df['start'] = headers['cds_start']
    cds_only_df['end'] = headers['cds_end']
    cds_only_df['id'] = headers['transcript_id']
    cds_only_df['score'] = 0
    cds_only_df['strand'] = headers['strand']

    bio.write_bed(cds_only_df, args.out)

if __name__ == '__main__':
    main()


