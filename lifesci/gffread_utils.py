""" This module includes parsing utilities for the gffread tool from cufflinks.
"""
import logging
logger = logging.getLogger(__name__)

import collections
import numpy as np
import pandas as pd
import re

import lifesci.fastx_utils as fastx_utils

# build up the header regular expression
transcript_id_re = r'(?P<transcript_id>\S+)'
gene_re = r"(?:gene=(?P<gene>\S+))?"
cds_re = r"(?:CDS=(?P<cds>\S+))?"
loc_re = r"loc:(?P<seqname>\S+)\|(?P<seq_pos>\S+)\|(?P<strand>[+-.])"
exons_re = r"exons:(?P<exons>\S*)"
segs_re = r"segs:(?P<segs>\S*)"

# there is no need for the initial ">" because pyfasta strips it
header_re = r"{}\s*{}\s*{}\s*{}\s*{}\s*{}".format(transcript_id_re,
                                                  gene_re,
                                                  cds_re,
                                                  loc_re,
                                                  exons_re,
                                                  segs_re)

header_re = re.compile(header_re)

default_novel_id_re = ""

fasta_header = collections.namedtuple(
        "fasta_header",
        "transcript_id,gene,cds_start,cds_end,seqname,strand,"
            "exon_starts,exon_ends,seg_starts,seg_ends"
)


def get_starts_ends(coords):
    """ This function parses gffread-style coordinate strings into a pair of arrays.
        The first array contains all of the start positions, and the second contains
        all of the ends.

        Args:
            coords (string) : a string of coordinate pairs. For example:
                1-56,57-261,262-543

        Returns:
            np.array : the first coordinate in each pair
            np.array : the second coordinate in each pair
    """

    s = re.split('-|,', coords)
    starts = np.array(s[0::2], dtype=int)
    ends = np.array(s[1::2], dtype=int)
    return starts, ends

def parse_header(header):
    """ This function parses out the exon (genomic) and segment (relative) coordinates
        from the gffread-style fasta header.

        Args:
            header (string) : a gffread-style header. For example:
                >ENST00000621500 gene=GPHB5 CDS=58-447 loc:14|63312835-63318879|- exons:63312835-63313116,63317646-63317850,63318824-63318879 segs:1-56,57-261,262-543

            exons_segs_re (compiled reg ex) : a compiled regular expression which
                parses the exons as the first "group" and the segments as the 
                second "group"

        Returns:
            fasta_header namedtuple, with the following fields:
                transcript_id (str, "ENST00000621500")
                gene (str, "GPHB5")
                cds_start (int, 58)
                cds_end (int, 447)
                seqname (str, "14")
                strand (str, "-")
                exon_starts (np.array of ints, the first (absolute) coordinate in each exon)
                exon_ends (np.array of ints, the second (absolute) coordinate in each exon)
                seg_starts (np.array of ints, the first (relative) coordinate in each exon)
                seg_ends (np.array of ints, the second (relative) coordinate in each exon)

    """
    global header_re

    m = header_re.match(header)

    if m is None:
        msg = "Failed parsing header. Header: '{}'".format(header)
        raise ValueError(msg)

    transcript_id = m.group('transcript_id')
    gene = m.group('gene')
    strand = m.group('strand')
    seqname = m.group('seqname')
    cds = m.group('cds')

    exons = m.group('exons')
    segs = m.group('segs')
    
    cds_start = None
    cds_end = None
    if cds is not None:
        cds_start, cds_end = get_starts_ends(cds)
        cds_start = cds_start[0]
        cds_end = cds_end[0]
    
    exon_starts, exon_ends = get_starts_ends(exons)
    seg_starts, seg_ends = get_starts_ends(segs)
    
    h = fasta_header(
        transcript_id=transcript_id,
        gene=gene,
        cds_start=cds_start,
        cds_end=cds_end,
        seqname=seqname,
        strand=strand,
        exon_starts=exon_starts,
        exon_ends=exon_ends,
        seg_starts=seg_starts,
        seg_ends=seg_ends
    )
    
    return h

def get_all_headers(transcripts_fasta_file):
    """ This function extracts the header information from all records in a
        given gffread-like fasta file. It returns them as a data frame. This
        could, for example, be used to extract only the CDS regions.

        Also, see the documentation for parse_header for more information.

        Args:
            transcripts_fasta_file (string or file-like): The path to the fasta file, or
                an open file-like object.
        Returns:
            pd.DataFrame: a data frame with the following fields:
                transcript_id (str, "ENST00000621500")
                gene (str, "GPHB5")
                cds_start (int, 58)
                cds_end (int, 447)
                seqname (str, "14")
                strand (str, "-")
                exon_starts (np.array of ints, the first (absolute) coordinate in each exon)
                exon_ends (np.array of ints, the second (absolute) coordinate in each exon)
                seg_starts (np.array of ints, the first (relative) coordinate in each exon)
                seg_ends (np.array of ints, the second (relative) coordinate in each exon)
    """
    msg = "Parsing all headers"
    logger.info(msg)

    transcripts = fastx_utils.get_read_iterator(transcripts_fasta_file)
    
    all_headers = []
    for transcript in transcripts:
        h = parse_header(transcript[0])
        all_headers.append(h)

    all_headers = [header for header in all_headers if header is not None]
    all_headers = pd.DataFrame(all_headers)
    all_headers = all_headers.fillna(value=0)

    all_headers['cds_start'] = all_headers['cds_start'].astype(int)
    all_headers['cds_end'] = all_headers['cds_end'].astype(int)

    return all_headers


def get_utr_info(fasta_record):
    """ Given a pyfasta entry from a gffread fasta file, this function extracts
        the 5' and 3' UTR lengths and sequences.

        Args:
            fasta_entry (tuple): a biopython-like fasta record, which has (at least)
                two elements. The first should be a gffread-like fasta header, and
                the second should be a DNA sequence.

                The sequence is not validated, but the header is parsed using the
                parse_header function in this package.

        Returns:
            None, if fasta_entry header does not include CDS information

                OR

            dictionary containing:
                transcript_id (string): the transcript identifer from the 
                    header (which is everything until the first space)

                five_prime_utr_len, three_prime_utr_len (ints): the length of
                    the respective UTRs

                five_prime_utr, three_prime_utr (strings): the respect UTR
                    sequences
    """
    header = parse_header(fasta_record[0])
    if header.cds_start is None:
        return None
        
    five_prime_utr_len = header.cds_start - 1
    five_prime_utr = fasta_record[1][:five_prime_utr_len]
    
    # we subtract three because gffread includes the stop codon
    three_prime_utr_len = header.seg_ends[-1] - header.cds_end - 3
    if three_prime_utr_len < 0:
        three_prime_utr_len = 0
        three_prime_utr = ""
    else:
        three_prime_utr = fasta_record[1][-1*three_prime_utr_len:]
    
    
    utr_info = {
        'transcript_id': header.transcript_id,
        'five_prime_utr_len': five_prime_utr_len,
        'three_prime_utr_len': three_prime_utr_len,
        'five_prime_utr': five_prime_utr,
        'three_prime_utr': three_prime_utr
    }
    
    return utr_info

def get_all_utrs(transcripts_fasta_file, progress_bar=False):
    """ This function extracts all of the annotated UTRs from the gffread-like
        transcripts file. In particular, UTRs are called as the parts of
        transcripts which occur before and after annotated CDSs. Transcripts
        with no annotated CDS region do not have UTRs according to this
        definition. The result is returned as a pandas data frame.

        Args:
            transcripts_fasta_file (string): the path to the gffread-like
                transcripts fasta file

            progress_bar (bool): whether to show a progress bar. Note that this
                requires iterating through the transcripts file twice to count
                the number of reads in the file.

        Returns:
            pd.DataFrame: a data frame with the following columns:
            
                transcript_id (string): the transcript identifer from the 
                    header (which is everything until the first space)

                five_prime_utr_len, three_prime_utr_len (ints): the length of
                    the respective UTRs

                five_prime_utr, three_prime_utr (strings): the respect UTR
                    sequences
    """
    if progress_bar:
        msg = "Counting sequences in file"
        logger.info(msg)
        read_count = fastx_utils.get_read_count(transcripts_fasta_file)
    else:
        read_count = None

    msg = "Extracting UTR sequences"
    logger.info(msg)

    transcripts = fastx_utils.get_read_iterator(transcripts_fasta_file)
    all_utr_info = parallel.apply_iter_simple(transcripts, get_utr_info, 
        progress_bar=progress_bar, total=read_count)

    all_utr_info = [utr_info for utr_info in all_utr_info if utr_info is not None]
    all_utr_info_df = pd.DataFrame(all_utr_info)

    return all_utr_info_df


