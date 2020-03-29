#! /usr/bin/env python3

import argparse

from Bio.Seq import translate

import lifesci.fastx_utils as fastx_utils

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script converts a DNA fasta file (i.e., with nucleic acids) "
        "into an amino acid fasta file (i.e., protein). It does not produce all reading "
        "frames, etc. It just converts the sequences as given.")

    parser.add_argument('dna', help="The DNA fasta file")
    parser.add_argument('out', help="The output (amino acid) file")

    parser.add_argument('--compress', help="If this flag is given, then the output will "
        "be compressed with gzip", action='store_true')
    
    args = parser.parse_args()

    records = fastx_utils.get_read_iterator(args.dna)
    protein_records = [
        (r[0], translate(r[1])) for r in records
    ]

    fastx_utils.write_fasta(protein_records, args.out, compress=args.compress)

if __name__ == '__main__':
    main()
