###
#   This module serves as a replacement for most pybedutils operations.
###

import collections
import csv
import logging
import os
import re
import shutil
import sys

import numpy as np
import pandas as pd
import tqdm

import pysam
from Bio.Seq import reverse_complement

import pyllars.collection_utils as collection_utils
import pyllars.math_utils as math_utils
import pyllars.pandas_utils as pd_utils
import pyllars.validation_utils as validation_utils
import pyllars.utils    


import lifesci.bam_utils as bam_utils
import lifesci.bio as bio
import lifesci.fastx_utils as fastx_utils

logger = logging.getLogger(__name__)

bed12_field_names = [
    "seqname", "start", "end", "id", "score", "strand",
    "thick_start", "thick_end", "color", "num_exons", "exon_lengths", "exon_genomic_relative_starts"
]

bed6_field_names = bed12_field_names[:6]
bed3_field_names= bed12_field_names[:3]

###
#   The following functions relate to bed file io. They are mostly wrappers
#   around pandas io.
###

def read_bed(filename, sep='\t', comment=None, header=None, 
                use_default_field_names=False, **kwargs):
    """ This function reads a bed file into a pandas data frame. By default, it 
        assumes the first line of the bed file actually gives the field names, 
        but they begin with hash marks (#). For example:

        #chr    #start  #end
        1   0   10

        Alternatively, the "default" bed field names can be used by setting the
        flag to True. Finally, a comment character (such as "#") and the field 
        names can be explicitly passed, in which case they will be sent to pandas.

        Args:
            filename (string): the path to the bed file

            header, comment, sep (strings): these arguments are passed through
                to pandas.read_csv.

            use_default_field_names (bool): If this flag is given, then the
                comment character will be used with the default bed field names.

            kwargs : these are passed through to the read_csv function


        Returns:
            pandas.DataFrame: a data frame with the bed information

        Imports:
            pandas
    """

    
    if header is None and not use_default_field_names:
        bed = pd.read_csv(filename, sep=sep, **kwargs)
        bed.columns = [c.replace("#", "") for c in bed.columns]

    else:
        if use_default_field_names:
            header = bed12_field_names

        bed = pd.read_csv(filename, sep=sep, header=None, comment=comment, **kwargs)
        num_columns = len(bed.columns)
        bed.columns = header[:num_columns]

    # either way, make sure the first column (seqname) is treated as a string
    seqname_column_name = bed.columns[0]
    seqname_column = bed[seqname_column_name]
    bed[seqname_column_name] = seqname_column.astype(str)

    return bed

def write_bed(data_frame, filename, compress=True, **kwargs):
    """ This function formats the data frame such that it will behave as
        a bed file. In particular, it writes a tab-delimited file and prepends
        the hash mark (#) to the field names in the header so they are treated
        as comments by typical bed software (such as bedtools).

        Args:
            data_frame (pandas.DataFrame) : a data frame representing bed objects

            filename (string) : the name of the output file

            compress (bool) : whether to gzip the output

            kwargs : these are passed through to the write_df function

        Returns:
            None

        Imports:
            misc.utils
            csv
            gzip (indirectly)
    """

    do_not_compress = not compress

    header = ['#{}'.format(c) for c in data_frame.columns]
    pd_utils.write_df(data_frame, filename, index=False, sep='\t', 
        header=header, do_not_compress=do_not_compress, quoting=csv.QUOTE_NONE, **kwargs)

def get_bed_df(bed):
    """ If bed is a string, treat it as a filename and read the bed file.
    Otherwise, return the object unchanged.

    Parameters
    ----------
    bed: object
        Presumably, this is either a bed-style data frame or the path to a bed
        file.

    Returns
    -------
    bed_df: bed-style pd.DataFrame
        If bed is a string, this is the bed df from the file at that location.
        Otherwise, this is just the input again.

    Raises
    ------
    ValueError: if bed is neither a string nor a pd.DataFrame
    """

    if isinstance(bed, str):
        bed = read_bed(bed)

    if not isinstance(bed, pd.DataFrame):
        msg = ("[get_bed_df]: bed does not appear to be a data frame. type: {}"
            .format(type(bed)))
        raise ValueError(msg)

    return bed

###
#   The following functions filter entries from bed files in one way or another.
###
def get_longest_features_by_end(bed_df):
    """ This function groups features by their strand-specific end coordinate
        (i.e., "end" for forward strand features and "start" for reverse strand
        features). It then selects the longest feature for each end coordinate.
        Ties are broken arbitrarily.

        Intuitively, this is used to remove truncated versions of features. It
        is similar in spirit to "bedtools merge" and "bedtools cluster."

        TODO : handle alternate field names

        Input:
            bed_df (pd.DataFrame) : a data frame corresponding to a bed file. It
                must contain the fields:
                    seqname
                    start
                    end
                    id
                    strand

        Returns:
            pd.DataFrame : a data frame with the truncated features removed

        Imports:
            pandas
    """
    
    # handle the separate strands
    m_forward = bed_df['strand'] == '+'
    m_reverse = bed_df['strand'] == '-'

    forward_orfs = bed_df[m_forward]
    reverse_orfs = bed_df[m_reverse]

    # group by the stop codon position
    forward_groups = forward_orfs.groupby(['seqname', 'end'])
    reverse_groups = reverse_orfs.groupby(['seqname', 'start'])

    # find the "start" position which results in the longest ORF
    forward_min_starts = forward_groups['start'].transform(min)
    reverse_max_ends = reverse_groups['end'].transform(max)

    # select only those ORFs
    m_longest_orf_in_forward_group = forward_orfs['start'] == forward_min_starts
    m_longest_orf_in_reverse_group = reverse_orfs['end'] == reverse_max_ends

    longest_orf_in_forward_group = forward_orfs[m_longest_orf_in_forward_group]
    longest_orf_in_reverse_group = reverse_orfs[m_longest_orf_in_reverse_group]

    # join the two strands back together and sort
    longest_orfs = pd.concat([longest_orf_in_forward_group, longest_orf_in_reverse_group])
    longest_orfs = longest_orfs.sort_values(['seqname', 'start'])
    longest_orfs = longest_orfs.reset_index(drop=True)
    
    # now, regroup the longest orfs to make sure we select each ORF only once
    id_groups = longest_orfs.groupby('id')
    idx_first = id_groups.apply(lambda x: x['start'].idxmax())
    
    longest_unique_orfs = longest_orfs.iloc[idx_first]
    
    return longest_unique_orfs

###
#   The following functions translate between transcript and genomic coordinates.
###

def parse_exon_start_end_length(exon, delimiter='-'):
    """ This function parses the start, end and length from an exon.
        It performs quite a bit of checking to ensure valid values.
        
        Example:
            3216021-3216967
        returns:
            (3216021,3216967,947)
            
        Args:
            exon (string) : the delimited start and stop of the exon
            
            delimiter (string) : the delimiter between the start and stop
            
        Returns:
            int : start
            int : end
            int : length
            
        Imports:
            misc.utils
    """         
    exon_start_end = exon.split(delimiter)
    if len(exon_start_end) != 2:
        msg = ("There was a problem parsing the exon coordinates: '{}'. Inner "
            "separator: '{}'".format(exon, delimiter))
        raise ValueError(msg)
        
    start = exon_start_end[0].strip()
    if not utils.is_int(start):
        msg = ("There was a problem parsing the exon coordinates: '{}'. Could "
            "not interpret \"start\" as an integer: '{}'".format(exon, exon_start_end[0]))
        raise ValueError(msg)
        
    end = exon_start_end[1].strip()
    if not utils.is_int(end):
        msg = ("There was a problem parsing the exon coordinates: '{}'. Could not "
            "interpret \"end\" as an integer: '{}'".format(exon, exon_start_end[0]))
        raise ValueError(msg)
        
    start = int(start)
    end = int(end)
    
    length = end-start+1
    if length < 0:
        msg = ("The length for the exon was negative: '{}'. Please ensure the coordinates "
                   "follow the BED conventions in which the \"smaller\" coordinate always "
                   "comes first, regardless of the strand of the feature.".format(exon))
        
        raise ValueError(msg)
        
    return (start, end, length)

def get_block_and_offset(pos, block_starts, block_lengths):
    """ This helper function finds the block and the offset within that block
        for the given position. This function verifies that the position falls
        within the selected block, based on length.

        This function is valid for both transcript (relative) positions and
        genomic positions; however, both "pos" and "block_starts" must be
        in the same coordinate system. That is, if "pos" is a relative
        position, then block_starts must also be in relative coordinates,
        and likewise with genomic coordinates.

        Args:
            pos (int): the position
            
            block_starts (np.array of ints): the start position of each block.
                These should follow bed conventions (i.e., monotonically 
                increasing regardless of strand)

            block_lengths (np.array of ints): the length of each block

        Returns:
            int: the block index
            int: the offset within that block

        Raises:
            ValueError if pos does not actually fall within one of the blocks

        Imports:
            numpy
    """
    # find the block of gen_pos
    # add one in case gen_pos is exactly on the start of a block
    i = np.searchsorted(block_starts, pos+1) - 1

    # find the offset within that block
    offset = pos - block_starts[i]

    # make sure it is valid
    if (offset < 0) or (offset > block_lengths[i]):
        msg = ("The pos does not fall within the block structure. pos: {}, "
            "block_starts: {}, block_lengths: {}".format(pos, block_starts,
            block_lengths))
        raise IndexError(msg)

    return i, offset


def get_relative_position(gen_pos, block_starts, block_lengths):
    """ This position calculates the relative position of a genomic position
        within a "block" structure (e.g., exons). That is, the function determines
        the block to which the genomic position belongs, and then the offset
        within that block. The relative position is the sum of the lengths of
        all preceding blocks plus the offset.

        The function also verifies that gen_pos actually falls within the block.
        So the offset within the block must be less than the length of the block.

        N.B. the returned value is 0-based wrt the blocks, and the block are
        considered BED-style half-open.

        Examples:
            block_starts = 10, 30, 50
            block_lengths = 10, 10, 10

            gen_pos = 10
            rel_pos = 0
            
            gen_pos = 35
            rel_pos = 15

            gen_pos = 50
            rel_pos = 20

        Args:
            gen_pos (int): the genomic position

            block_starts (np.array): the start position of each block. These
                *must* be sorted in increasing order.

            block_lengths (np.array): the length of the blocks, where the order
                corresponds to that of block_starts

        Returns:
            rel_pos (int): the relative position within the block structure

        Imports:
            numpy

        Raises:
            IndexError if gen_pos does not actually fall within one of the blocks
    """
    i, offset = get_block_and_offset(gen_pos, block_starts, block_lengths)
    rel_pos = np.sum(block_lengths[:i]) + offset

    return rel_pos

def get_gen_pos(rel_pos, start, block_lengths, block_starts, block_relative_starts):
    """ This function converts from a relative position within a bed block
        structure to its genomic position. Roughly, this function determines
        the block to which the rel_pos belongs and then determines the 
        genomic position of the start of that block. The offset within the
        block then gives the genomic position of rel_pos.
        
        N.B. The blocks are considered bed-style half-open.
        
        Examples:
            start = 10
            block_lengths = 5,10,5,10,5
            block_starts = 0,5,15,20,30
            block_relative_starts = 0,15,35,50,80
            
            rel_pos: 2
            gen_pos: 12
            
            rel_pos: 16
            gen_pos: 46
            
        Args:
            rel_pos (int): the position within a block structure (i.e., transcript)
            
            start (int): the genomic start of the transcript
            
            block_lengths (np.array of ints): the length of each block
            
            block_starts (np.array of ints): the start of each block within the
                transcript structure. This is probably calculated somehow like:
                
                    block_starts = np.zeros(len(block_lengths), dtype=int)
                    block_starts[1:] = np.cumsum(block_lengths)[:-1]
            
            block_relative_starts (np.array of ints): the (genomic) relative 
                block starts (12th bed column)
            
        Returns:
            int: the genomic position
            
        Raises:
            IndexError if the rel_pos is negative or after the end of the last block
            
        Imports:
            numpy
            misc.math_utils
    """    
    # first, determine the block into which this position falls, and the offset within that block
    block, offset = get_block_and_offset(rel_pos, block_starts, block_lengths)
    
    # this is too much debug information; unit tests are designed to test this
    #msg = "offset: {}".format(offset)
    #logger.debug(msg)
    
    # get the genomic start of the block
    block_gen_start = start + block_relative_starts[block]
    
    #msg = "block_gen_start: {}".format(block_gen_start)
    #logger.debug(msg)
    
    # and shift in the block
    gen_pos = block_gen_start + offset

    #msg = "gen_pos: {}".format(gen_pos)
    
    return gen_pos

def convert_genomic_coords_to_bed_blocks(coords, outer_sep=',', inner_sep='-'):
    """ This function converts a list of exon coordinates into the "blocks"
        columns for the a BED12 file ("num_exons", "exon_lengths", 
        "exon_genomic_relative_starts").
        
        Example:
            3216021-3216967, 3421701-3421900, 3670551-3671347
        becomes:
            {
                num_exons: 3,
                exon_lengths: "947,200,797"
                exon_genomic_relative_starts: "0,205680,454530"
            }
            
        N.B. Because the BED block coordinates are relative to the start of the
                first exon, base-0 vs. base-1 is not a problem.
                
        N.B. The coordinates are taken to be INCLUSIVE (e.g., as with GTF).
                
        Args:
            coords (string) : the string from which the coordinates will be parsed
            
            outer_sep (string) : the separator for each exon ("," in the example)
            
            inner_sep (string) : the separator for the start and end of each exon 
                ("-" in the example)
            
        Returns:
            dict : a dictionary containing:
                num_exons
                exon_lengths
                exon_genomic_relative_starts
                
            Presumably, these would be used to create columns in a data frame or BED12 file.
            
        Imports:
            misc.utils (indirectly)
            misc.math_utils
    """
    exons = coords.split(outer_sep)
    
    num_exons = len(exons)
    exon_lengths = []
    relative_starts = []
    
    if num_exons < 1:
        msg = ("There was a problem parsing the outer coordinates: '{}'. Outer "
            "separator: '{}'".format(coords, outer_sep))
        raise ValueError(msg)
        
    # we need to handle the first exon specially because we need its start
    (start, end, length) = parse_exon_start_end_length(exons[0], inner_sep)
    
    first_start = start
    exon_lengths.append(length)
    relative_starts.append(0)
    
    for exon in exons[1:]:
        (start, end, length) = parse_exon_start_end_length(exon, inner_sep)
        
        relative_start = start - first_start
        exon_lengths.append(length)
        relative_starts.append(relative_start)
        
    # now, make sure the relative starts are strictly increasing
    if not math_utils.is_monotonic(relative_starts):
        msg = ("There was a problem parsing the coordinates: '{}'. The outer "
            "coordinates are not monotonically increasing.".format(coords))
        raise ValueError(msg)
        
    # we need the values as csv
    exon_lengths = ','.join(str(l) for l in exon_lengths)
    relative_starts = ','.join(str(s) for s in relative_starts)
        
    return {
        'num_exons': num_exons,
        'exon_lengths': exon_lengths,
        'exon_genomic_relative_starts': relative_starts
    }

###
#   The following functions are utilities for working with bed12 files.
###

def get_bed_12_feature_length(bed_feature):
    """ This function calculates the length of a split BED12 feature by summing
        the lengths of each exon.

        Input:
            bed_feature (pd.Series) : a BED12 feature. It must contain the field:
                exon_lengths

        Returns:
            int : the length of all exons
    """
    lengths = str(bed_feature['exon_lengths']).split(",")
    total_length = sum(int(l) for l in lengths)
    return total_length

def split_bed12_blocks(bed_entry):
    """ Given a BED12+ entry using the bio.bed12_field_names, this function
        splits the entry into separate entries for each of the defined blocks.
        It returns a new BED12+ entry as a dictionary.

        Args:
            bed_entry (dict-like): a BED12+ entry which can be indexed as a
                dictionary (e.g., a dict or pd.Series). The following field 
                names must match those in bio.bed12_field_names:

                    seqname, start, id, score, strand, exon_lengths,
                    exon_genomic_relative_starts

        Returns:
            list of dicts: a list of bed6+2 entries of the blocks. The
                additional fields are:

                    exon_index: the 0-based index of the exon. This is always
                        relative to "start". So, for reverse-strand features,
                        exon 0 is actually the "last" exon in the feature

                    transcript_start: the 0-based index of the start of this
                        exon within the transcript. As with exon_index, this
                        is relative to "start". 

        Imports:
            numpy
    """
    exon_lengths = np.fromstring(bed_entry['exon_lengths'], 
        sep=",", dtype='int')
    exon_rel_starts = np.fromstring(bed_entry['exon_genomic_relative_starts'], 
        sep=",", dtype='int')

    exon_starts = bed_entry['start'] + exon_rel_starts
    exon_ends = exon_starts + exon_lengths

    transcript_starts = np.zeros(len(exon_starts)+1, dtype=int)
    transcript_starts[1:] = np.cumsum(exon_lengths)

    exons = []

    for i, (start, end, transcript_start) in enumerate(zip(exon_starts, exon_ends, transcript_starts)):
        exon = {
            "seqname": bed_entry["seqname"],
            "start": start,
            "end": end,
            "id": bed_entry["id"],
            "score": bed_entry["score"],
            "strand": bed_entry["strand"],
            "exon_index": i,
            "transcript_start": transcript_start
        }
        
        exons.append(exon)

    return exons

def split_bed12(bed_df, num_cpus=1, progress_bar=False):
    """ This function splits a BED12+ data set, represented as a pandas dataframe,
        into a BED6+ data set. Each block (i.e., exon) from the BED12+ data set
        is an individual feature in the result. The results are sorted.

        Args:
            bed_df (pd.DataFrame): a BED12+ data set. The following field 
                names must match those in bio.bed12_field_names:

                    seqname, start, id, score, strand, exon_lengths,
                    exon_genomic_relative_starts

            num_cpus (int): the number of CPUs to parallelize across

            progress_bar (bool): whether to use a tqdm progress bar

        Returns:
            bed6_df (pd.DataFrame): a BED6+2 dataset. The additional fields are:

                exon_index: the 0-based index of the exon. This is always
                    relative to "start". So, for reverse-strand features,
                    exon 0 is actually the "last" exon in the feature

                transcript_start: the 0-based index of the start of this
                    exon within the transcript. As with exon_index, this
                    is relative to "start". 
    """    
    columns = bed6_field_names + ['exon_index', 'transcript_start']

    if len(bed_df) == 0:
        msg = ("[bed_utils.split_bed12]: Attempting to split empty bed data "
            "frame")
        logger.warning(msg)
        return pd.DataFrame(columns=columns)
    
    exons = parallel.apply_parallel(
        bed_df,
        num_cpus,
        split_bed12_blocks,
        progress_bar=progress_bar)

    exons_flat = collection_utils.flatten_lists(exons)
    exons_df = pd.DataFrame(exons_flat)


    return exons_df[columns]

###
#   The following functions are used to trim away various parts bed12 entries.
###

def retain_thick_only(bed_entry, inplace=False):
    """ Given a BED12+ entry using the bio.bed12_field_names, this function
        trims away all of the entry except for the "thick" part. It returns
        a new BED12+ entry as a dictionary.
        
        This function is not implemented particularly efficiently; however, it
        is implemented reasonably clearly in an attempt to minimize off-by-one
        and similar errors.

        Optionally, the update can happen inplace.
        
        Args:
            bed_entry (dict-like): a BED12+ entry which can be indexed as a
                dictionary (e.g., a dict or pd.Series). The following field 
                names must match those in bio.bed12_field_names:
                
                    start, end, thick_start, end, thick_end
                    num_exons, exon_lengths, exon_genomic_relative_starts
                    
                Additionally, the object must have a "copy" function.

            inplace (bool): whether to create a new version of the object or
                update the existing one. If this is True, then the object does
                not need a "copy" function.
                    
        Returns:
            dict-like: an object of the same type as bed_entry with the fields 
                updated to include only the "thick" part of the original 
                record. All other fields are passed through unchanged. The
                updated version is returned regardless of the value of inplace.
                
        Imports:
            numpy
    """
    # first, extract the relevant values
    rel_starts = bed_entry['exon_genomic_relative_starts']
    rel_starts = np.array(rel_starts.split(","), dtype=int)
    
    lengths = bed_entry['exon_lengths']
    lengths = np.array(lengths.split(","), dtype=int)
    
    start = bed_entry['start']
    end = bed_entry['end']
    
    thick_start = bed_entry['thick_start']
    thick_end = bed_entry['thick_end']

    # find genomic starts
    genomic_starts = start + rel_starts

    # find exon of thick_start
    # add one to thick_start in case it is exactly on an exon start
    i = np.searchsorted(genomic_starts, thick_start+1) - 1

    # find offset within that exon
    offset = thick_start - genomic_starts[i]

    # update the length of that exon
    lengths[i] -= offset

    # discard the previous lengths
    lengths = lengths[i:]

    # discard the previous relative starts
    genomic_starts = genomic_starts[i:]

    # update the first genomic start to point to the thick_start
    genomic_starts[0] = thick_start

    # transform the genomic starts to relative starts by subtracting thick_start
    genomic_starts = genomic_starts - thick_start

    rel_starts = genomic_starts
    start = thick_start

    # now, remove everything after thick_end

    # find the genomic_starts
    genomic_starts = rel_starts + start

    # find the exon of thick_end
    i = np.searchsorted(genomic_starts, thick_end) - 1

    # find the offset within that exon
    offset = thick_end - genomic_starts[i]

    # update the length of that exon
    lengths[i] = offset

    # discard the remaining exons
    lengths = lengths[:i+1]
    genomic_starts = genomic_starts[:i+1]

    # transform back to relative coordinates
    rel_starts = genomic_starts - start
    
    # and update the number of exons
    num_exons = len(rel_starts)
    
    if inplace:
        ret = bed_entry
    else:
        ret = bed_entry.copy()
    
    ret['start'] = thick_start
    ret['end'] = thick_end
    
    ret['num_exons'] = num_exons
    ret['exon_lengths'] = ",".join(str(l) for l in lengths)
    ret['exon_genomic_relative_starts'] = ",".join(str(s) for s in rel_starts)
    
    return ret

def get_no_thick_mask(bed_df, no_thick_value=-1):
    """ This helper function creates a mask for the bed12+ data frame that
        indicates all entries which do not have "thick" coordinates.

        Args:
            bed_df (pd.DataFrame): the bed12+ data frame

            no_thick_value (obj): the value used to check that the thick value
                is missing

        Returns:
            bool mask: a mask in which "True" indicates the entry at
                that position has no "thick" coordinates
    """
    m_no_thick_start = bed_df['thick_start'] == -1
    m_no_thick_end = bed_df['thick_end'] == -1
    m_no_thick = m_no_thick_start & m_no_thick_end
    return m_no_thick

def retain_all_no_thick(bed_df, no_thick_value=-1):
    """ Given a bed12+ data frame, this function removes all entries which have
        a "thick" part defined. It returns a copy with only the transcripts with
        no thick part.

        Args:
            bed_df (pd.DataFrame): the bed12+ data frame

            no_thick_value (obj): the value used to check that the thick value
                is missing

        Returns:
            pd.DataFrame: a copy of bed_df which includes only those entries
                which did not have a thick part defined            
    """
    m_no_thick = get_no_thick_mask(bed_df, no_thick_value=no_thick_value)
    ret = bed_df[m_no_thick]
    ret = ret.copy()
    return ret

def retain_all_thick_only(bed_df, prune_no_thick=True, progress_bar=True, num_cpus=1):
    """ Given a bed12+ data frame, this function removes the "non-thick" part
        from all of the entries. Optionally, it will also prune entries which
        do not have a thick part. This is detected by "-1" values in both
        "thick" fields.

        This function is largely a wrapper aroiund retain_thick_only, so see
        the documentation of that function for more details.

        This function always returns a copy.

        Args:
            bed_df (pd.DataFrame): a data frame which uses bio.bed12_field_names

            prune_no_thick (bool): whether to remove entries with no thick part.
                If False, they will still be in the data frame, but all
                coordinates will be "-1".

            progress_bar (bool): whether to show a progress bar

            num_cpus (int): the number of CPUs to use

        Returns:
            pd.DataFrame: a copy of the original data frame with the bed entries
                updated to include only the thick parts. Depending on the value
                of prune_no_thick, the entries with no thick part will be
                removed.
    """
    if prune_no_thick:
        m_no_thick = get_no_thick_mask(bed_df)
    else:
        m_no_thick = np.zeros(len(bed_df), dtype=bool)

    ret = parallel.apply_parallel(bed_df[~m_no_thick], 
                            num_cpus, 
                            retain_thick_only, 
                            progress_bar=progress_bar)

    ret = pd.DataFrame(ret)
    return ret


def retain_all_five_prime_of_thick(bed_df:pd.DataFrame, progress_bar:bool=True,
        num_cpus:int=1):
    """ Given a bed12+ data frame, this function retains only the portion of
    each entry that is 5' of the annotated "thick" region. It also removes
    all entries which do not have a "thick" region.

    N.B. This function *does* consider the strand.

    This function removes all 0-length 5' leaders, as well. However, it can
    still happen that an exon in the 5' leader has length 0 if the "thick"
    portion begins exactly at an exon start.

    This function is largely a wrapper around retain_all_thick_only, so
    please see that documentation for more details.

    This function always returns a copy.

    Parameters
    ----------
    bed_df: pd.DataFrame
        a data frame which uses bio.bed12_field_names

    progress_bar: bool
        whether to show a progress bar

    num_cpus: int
        the number of CPUs to use

    Returns
    -------
    five_prime_only: pd.DataFrame
        a copy of the original data frame with the bed entries updated to 
        include only the regions 5' of the thick parts.
    """
    # first, make sure we remove everything without a thick region

    # we will be updating these values, so it is important we keep
    # track of these first
    m_no_thick = get_no_thick_mask(bed_df)

    # also, get those which do not have a 5' leader
    m_forward = bed_df['strand'] == '+'
    m_reverse = bed_df['strand'] == '-'

    m_equal_starts = bed_df['start'] == bed_df['thick_start']
    m_equal_ends = bed_df['end'] == bed_df['thick_end']

    m_no_leader = (m_forward & m_equal_starts) | (m_reverse & m_equal_ends)

    m_ignore = m_no_thick | m_no_leader

    bed_df = bed_df[~m_ignore].copy()

    m_forward = bed_df['strand'] == '+'
    m_reverse = bed_df['strand'] == '-'

    bed_df.loc[m_forward, 'thick_end'] = bed_df.loc[m_forward, 'thick_start']
    bed_df.loc[m_forward, 'thick_start'] = bed_df.loc[m_forward, 'start']
    
    bed_df.loc[m_reverse, 'thick_start'] = bed_df.loc[m_reverse, 'thick_end']
    bed_df.loc[m_reverse, 'thick_end'] = bed_df.loc[m_reverse, 'end']

    ret = retain_all_thick_only(bed_df, progress_bar=progress_bar, 
            num_cpus=num_cpus)

    if len(ret) == 0:
        msg = ("[bed_utils.retain_all_five_prime_of_thick]: No annotated 5' "
            "regions before the \"thick\" segments")
        logger.warning(msg)
        ret = pd.DataFrame(columns=bed_df.columns.copy())

    return ret

def retain_all_three_prime_of_thick(bed_df:pd.DataFrame, 
        progress_bar:bool=True, num_cpus:int=1):
    """ Given a bed12+ data frame, this function retains only the portion of
    each entry that is 3' of the annotated "thick" region. It also removes
    all entries which do not have a "thick" region.

    N.B. This function *does* consider the strand.

    This function is largely a wrapper around retain_all_thick_only, so
    please see that documentation for more details.

    This function always returns a copy.

    Parameters
    ----------
    bed_df: pd.DataFrame
        a data frame which uses bio.bed12_field_names

    progress_bar: bool
        whether to show a progress bar

    num_cpus: int
        the number of CPUs to use

    Returns
    -------
    three_prime_only: pd.DataFrame
        a copy of the original data frame with the bed entries updated to 
        include only the regions 3' of the thick parts.
    """
    # first, make sure we remove everything without a thick region

    # we will be updating these values, so it is important we keep
    # track of these first
    m_no_thick = get_no_thick_mask(bed_df)

    # also, get those which do not have a 3' leader
    m_forward = bed_df['strand'] == '+'
    m_reverse = bed_df['strand'] == '-'

    m_equal_starts = bed_df['start'] == bed_df['thick_start']
    m_equal_ends = bed_df['end'] == bed_df['thick_end']

    m_no_trailer = (m_forward & m_equal_ends) | (m_reverse & m_equal_starts)

    m_ignore = m_no_thick | m_no_trailer

    bed_df = bed_df[~m_ignore].copy()

    m_forward = bed_df['strand'] == '+'
    m_reverse = bed_df['strand'] == '-'

    bed_df.loc[m_forward, 'thick_start'] = bed_df.loc[m_forward, 'thick_end']
    bed_df.loc[m_forward, 'thick_end'] = bed_df.loc[m_forward, 'end']
    
    bed_df.loc[m_reverse, 'thick_end'] = bed_df.loc[m_reverse, 'thick_start']
    bed_df.loc[m_reverse, 'thick_start'] = bed_df.loc[m_reverse, 'start']

    ret = retain_all_thick_only(bed_df, progress_bar=progress_bar, 
            num_cpus=num_cpus)

    if len(ret) == 0:
        msg = ("[bed_utils.retain_all_three_prime_of_thick]: No annotated 3' "
            "regions after the \"thick\" segments")
        logger.warning(msg)
        ret = pd.DataFrame(columns=bed_df.columns.copy())

    return ret


def retain_before_thick_only(bed_entry, inplace=False):
    """ Given a bed12+ entry using the bio.bed12_field_names, this function
        trims away all of the entry except for the bit before the "thick" part. 
        It returns a new bed12+ entry as a dictionary-like.
        
        Because the updated bed entry does not contain any of the original 
        "thick" part, "thick_start" and "thick_end" are set to -1.
        
        This function essentially wraps retain_thick_only, so its efficiency is
        similar to that of retain_thick_only (i.e., not so efficient).

        Optionally, the update can happen inplace.
        
        N.B. "before" in this context always refers to "smaller genomic
        coordinates." Thus, for forward strand entries, this refers to the
        region 5' of the "thick" part, while for negative strand entries,
        this refers to the 3' region relative to the "thick" part.
        
        Args:
            bed_entry (dict-like): a bed12+ entry which can be indexed as a
                dictionary (e.g., a dict or pd.Series). The following field 
                names must match those in bio.bed12_field_names:
                
                    start, end, thick_start, end, thick_end
                    num_exons, exon_lengths, exon_genomic_relative_starts
                    
                Additionally, the object must have a "copy" function.

            inplace (bool): whether to create a new version of the object or
                update the existing one. If this is True, then the object does
                not need a "copy" function.
                    
        Returns:
            dict-like: an object of the same type as bed_entry with the fields 
                updated to include only the bit before the "thick" part of the 
                original entry. All other fields are passed through unchanged.
                The updated version is returned regardless of the value of 
                inplace.
    """
    if inplace:
        ret = bed_entry
    else:
        ret = bed_entry.copy()
    
    new_thick_start = bed_entry['start']
    new_thick_end = bed_entry['thick_start']
    
    ret['thick_start'] = new_thick_start
    ret['thick_end'] = new_thick_end
    
    retain_thick_only(ret, inplace=True)
    
    ret['thick_start'] = -1
    ret['thick_end'] = -1
    
    return ret


def retain_after_thick_only(bed_entry, inplace=False):
    """ Given a bed12+ entry using the bio.bed12_field_names, this function
        trims away all of the entry except for the bit after the "thick" part. 
        It returns a new bed12+ entry as a dictionary-like.
        
        Because the updated bed entry does not contain any of the original 
        "thick" part, "thick_start" and "thick_end" are set to -1.
        
        This function essentially wraps retain_thick_only, so its efficiency is
        similar to that of retain_thick_only (i.e., not so efficient).

        Optionally, the update can happen inplace.
        
        N.B. "after" in this context always refers to "larger genomic
        coordinates." Thus, for forward strand entries, this refers to the
        region 3' of the "thick" part, while for negative strand entries,
        this refers to the 5' region relative to the "thick" part.
        
        Args:
            bed_entry (dict-like): a bed12+ entry which can be indexed as a
                dictionary (e.g., a dict or pd.Series). The following field 
                names must match those in bio.bed12_field_names:
                
                    start, end, thick_start, end, thick_end
                    num_exons, exon_lengths, exon_genomic_relative_starts
                    
                Additionally, the object must have a "copy" function.

            inplace (bool): whether to create a new version of the object or
                update the existing one. If this is True, then the object does
                not need a "copy" function.
                    
        Returns:
            dict-like: an object of the same type as bed_entry with the fields 
                updated to include only the bit after the "thick" part of the 
                original entry. All other fields are passed through unchanged.
                The updated version is returned regardless of the value of 
                inplace.
    """

    if inplace:
        ret = bed_entry
    else:
        ret = bed_entry.copy()
    
    new_thick_start = bed_entry['thick_end']
    new_thick_end = bed_entry['end']
    
    ret['thick_start'] = new_thick_start
    ret['thick_end'] = new_thick_end
    
    retain_thick_only(ret, inplace=True)
    
    ret['thick_start'] = -1
    ret['thick_end'] = -1
    
    return ret

###
#   The following functions are for merging bed entries.
###

def merge_intervals(interval_starts, interval_ends, interval_info=None):
    """ This function merges a given list of intervals based on overlaps. It
        merges intervals which overlap by at least 1bp. This function follows
        BED semantics, so the starts are inclusive and the ends are exclusive.

        Internally, a variant of the "chrom-sweep" algorithm is used to
        efficiently find the intersections. See [Pedersen et al., Genome Biology
        17:118, 2016] for a good description of this algorithm.

        Optionally, arbitrary information can be attached with each interval. If
        this information is provided, then the information for all intervals which
        are merged will be included as a list.

        The function returns a list of all merged interval start and ends.

        Algorithmically, this function is efficient, O(n), where n is the number
        of intervals. However, it does create some intermediate copies of the 
        data while, for example, sorting.

        Args:
            interval_starts, interval_ends (list-likes): a list of the start
                and end positions of the intervals. The arrays should match,
                so interval_starts[i] matches to interval_ends[i].

            interval_info (list-like, or None): a list of additional information
                about the intervals, or None if only the positions of the merged
                intervals are desired.

            N.B. The intervals will be sorted as the first step of this function,
            so they do not need to be pre-sorted.

        Returns:
            np.array: starts
            np.array: ends
            np.array of lists: interval_info for the merged intervals

            The sorted starts and ends of the merged intervals, as well as a list
            of the information about the merged intervals, if given.
    """
    if len(interval_starts) == 0:
        return [], [], []

    interval_starts = np.array(interval_starts)
    interval_ends = np.array(interval_ends)

    sorted_interval_start_indices = np.argsort(interval_starts)
    
    interval_starts = interval_starts[sorted_interval_start_indices]
    interval_ends = interval_ends[sorted_interval_start_indices]

    if interval_info is not None:
        interval_info = np.array(interval_info)
        interval_info = interval_info[sorted_interval_start_indices]
        next_interval_info = interval_info[0]

    else:
        next_interval_info = None

    next_interval = 0
    num_intervals = len(interval_starts)
    next_interval_start = interval_starts[0]
    next_interval_end = interval_ends[0]

    cache = []

    merged_starts = []
    merged_ends = []
    merged_info = []

    cur_interval_start = None
    cur_interval_end = None
    cur_interval_info = None

    while next_interval < num_intervals:
        # first, remove everything from the cache that ended before this begins
        cache = [c for c in cache if interval_ends[c] > next_interval_start]
        
        # if we cleared the cache, then we just finished an interval
        if len(cache) == 0:
            if cur_interval_start is not None:
                # add it to the output list
                merged_starts.append(cur_interval_start)
                merged_ends.append(cur_interval_end)
                merged_info.append(cur_interval_info)
                
            # and start a new interval
            cur_interval_start = next_interval_start
            cur_interval_end = next_interval_end
            cur_interval_info = [next_interval_info]
            
        else:
            # otherwise, extend the previous interval
            if next_interval_end > cur_interval_end:
                cur_interval_end = next_interval_end
                cur_interval_info.append(next_interval_info)
            
        # add the next interval to the current cache
        cache.append(next_interval)
        
        # and advance
        next_interval += 1
        if next_interval < num_intervals:

            next_interval_start = interval_starts[next_interval]
            next_interval_end = interval_ends[next_interval]

            if interval_info is not None:
                next_interval_info = interval_info[next_interval]
            
    # add the last merged interval to the output list
    merged_starts.append(cur_interval_start)
    merged_ends.append(cur_interval_end)
    merged_info.append(cur_interval_info)

    merged_starts = np.array(merged_starts)
    merged_ends = np.array(merged_ends)
    return merged_starts, merged_ends, merged_info

def merge_all_intervals(bed, split=False):
    """ This function merges all intervals in the given BED6+ data frame. It
        returns the merged intervals and the "id" of the intervals which
        overlap.

        N.B. By default, this function *does not* split BED12+ records. It 
            asssumes they are already split, and the "id" field gives useful 
            information. If specified, though, the function *will* split them
            first.

        This function is largely a wrapper around merge_intervals, so please
        see the documentation for that function for more information.

        Args:
            bed (pd.DataFrame): a data frame corresponding to a BED6+ data set

        Returns:
            pd.DataFrame: a BED6+1 data frame with the following fields:
                seqname (str): the seqname of the merged interval
                
                start (int): the (closed) start of the merged interval

                end (int): the (open) end of the merged interval

                id (string): the IDs of the merged intervals, separated by three
                    colons (":::")

                score (int): the number of original intervals in the merged interval
                
                strand: the strand of the merged interval

                merged_ids (list of string): a list of the original identifiers
                    of the intervals which make up the merged interval
    """
    if split:
        msg = "Splitting the BED12 records"
        logger.debug(msg)

        bed = split_bed12(bed)
    
    seqnames = bed['seqname'].unique()
    strands = ("+", "-")
    
    all_matches = []
    
    for strand in strands:
        m_strand = bed['strand'] == strand

        for seqname in seqnames:
            m_seqname = bed['seqname'] == seqname

            m_filter = m_strand & m_seqname

            if sum(m_filter) == 0:
                continue
        
            interval_starts = bed.loc[m_filter, 'start']
            interval_ends = bed.loc[m_filter, 'end']
            interval_info = bed.loc[m_filter, 'id']

            merged_intervals = merge_intervals(interval_starts, 
                                                interval_ends, interval_info)

            (merged_starts, merged_ends, merged_info) = merged_intervals

            merged_ids_str = [":::".join(merged_id) for merged_id in merged_info]
            num_merged_ids = [len(merged_id) for merged_id in merged_info]

            
            df = pd.DataFrame()
            df['start'] = merged_starts
            df['end'] = merged_ends
            df['seqname'] = seqname
            df['strand'] = strand
            df['score'] = num_merged_ids
            df['id'] = merged_ids_str
            df['merged_ids'] = merged_info
            
            all_matches.append(df)
            
    all_matches = pd.concat(all_matches)
    
    fields = bed6_field_names + ['merged_ids']
    return all_matches[fields]

position_interval_intersection = collections.namedtuple(
    "position_interval_intersection",
    "interval_index,relative_offset,interval_info,position_info"
)

###
#   The following functions are for finding the intersections of points and
#   intervals based on genomic coordinates. The results include both the
#   objects which overlap and the relative position within the transcript.
###

def get_position_intersections(positions, interval_starts, interval_ends,
            interval_info=None, position_info=None):

    """ This function finds the intersections of a set of (1bp) points and a
        set of intervals, specified by (inclusive) start and (exclusive) end
        positions. Furthermore, it allows arbitrary information to be attached
        with each of the intervals and each of the positions.

        Internally, a variant of the "chrom-sweep" algorithm is used to
        efficiently find the intersections. See [Pedersen et al., Genome Biology
        17:118, 2016] for a good description of this algorithm.

        This function returns a list of all intersections, including the interval
        index, the relative position of the intersection within that index, and
        the arbitrary information associated with that interval and position.

        Algorithmically, this function is efficient. However, it does create some
        intermediate copies of the data while, for example, sorting.

        Args:
            p_site_positions (list-like): a list of the 1bp point positions

            interval_starts, interval_ends (list-likes): a list of the start
                and end positions of the intervals. The arrays should match,
                so interval_starts[i] matches to interval_ends[i].

            interval_info (list-like): a list of additional information associated
                with the intervals. The container must be index-able, and it should
                match with the interval_starts and interval_ends. So 
                interval_info[i] should return the information associated with
                interval_start[i].

            position_info (list-like): a list of additional information associated
                with the positions. The container must be index-able, and it should
                match with the positions. So position_info[i] should return the 
                information associated with
                position[i].

            N.B. The intervals will be sorted as the first step of this function,
            so they do not need to be pre-sorted.

        Returns:
            list of 4-tuples: Each intersection generates a four-tuple with
                the following pieces of data:

                index 0: the index of the interval

                index 1: the relative position of the intersection within the
                    interval. For example, if a p_site intersects the first base
                    within an interval, this value is 0.

                index 2: the additional information associated with this interval,
                    or None if no additional information was given.

                index 3: the additional information associated with this
                    position, or None if no additional information was given.
    """
    matches = []
    cache = []

    num_positions = len(positions)
    num_intervals = len(interval_starts)
    
    # TODO: this is wasteful
    if position_info is None:
        position_info = [None] * num_positions

    if interval_info is None:
        interval_info = [None] * num_intervals

    positions = np.array(positions, dtype=int)
    interval_starts = np.array(interval_starts, dtype=int)
    interval_ends = np.array(interval_ends, dtype=int)
    position_info = np.array(position_info, dtype=object)
    interval_info = np.array(interval_info, dtype=object)

    sorted_position_indices = np.argsort(positions)
    positions = positions[sorted_position_indices]
    position_info = position_info[sorted_position_indices]
    
    sorted_interval_start_indices = np.argsort(interval_starts)
    interval_starts = interval_starts[sorted_interval_start_indices]
    interval_ends = interval_ends[sorted_interval_start_indices]
    interval_info = interval_info[sorted_interval_start_indices]
    
    positions = np.append(positions, np.inf)
    interval_starts = np.append(interval_starts, np.inf)

    next_p_site_index = 0
    next_p_site_position = positions[0]

    next_exon_index = 0
    next_exon_start = interval_starts[0]

    while next_p_site_position != np.inf:

        # do we grab the p_site or the exon
        if next_p_site_position < next_exon_start:

            # then we take the p_site

            # first, remove everything from the cache which ends before this
            # technically, keep everything which ends after this

            # we use ">" rather than ">=" because our intervals are open on
            # the "end" side
            cache = [c for c in cache if interval_ends[c] > next_p_site_position]

            # now, get the relative positions of everything in the cache
            rel_pos = [position_interval_intersection(c, next_p_site_position - interval_starts[c], 
                        interval_info[c], position_info[next_p_site_index])
                          for c in cache]

            # and add them to the output
            if len(rel_pos) > 0:
                matches.append(rel_pos)

            # advance to the next p-site
            next_p_site_index += 1

            next_p_site_position = positions[next_p_site_index]

        else:
            # just add this exon to the cache
            cache.append(next_exon_index)

            # advance to the next exon
            next_exon_index += 1

            next_exon_start = interval_starts[next_exon_index]
            
    matches = collection_utils.flatten_lists(matches)
    #matches = np.array(matches, dtype=position_interval_intersection)
    return matches

def get_all_position_intersections(positions_bed, intervals_bed, logger=logger):
    """ This function finds all intersections of the positions in the intervals.
        It expects both inputs to be BED6(+) data frames. In particular:

        * IT TAKES THE START OF ENTRIES IN POSITIONS_BED AS THE POSITIONS FOR
          FINDING INTERSECTIONS *

        * IT DOES NOT SPLIT INTERVALS_BED BASED ON BED12 FIELDS *

        This function is largely a wrapper around get_position_intersections,
        so please see the documentation of that function for more details.

        Args:
            positions_bed (pd.DataFrame): a bed6+ data frame. N.B., the start
                field is taken as the position for finding intersections.

            intervals_bed (pd.DataFrame): a bed6+ data frame.

            logger (logging.logger): optionally, a non-default logger can be
                given, in which case it will be used for logging. Otherwise,
                the lifesci.bed_utils logger will be used.

        Returns:
            A list of named 4-tuples with the following names and values:

                interval_index (index 0): the index of the interval

                relative_offset (index 1): the relative position of the intersection 
                    within the interval. For example, if a p_site intersects the 
                    first base within an interval, this value is 0.

                interval_info (index 2): the bed entry for the matching interval
                
                position_info (index 3): the bed entry for the matching position

            The type of the tuples is: bed_utils.position_interval_intersection
    """
    seqnames = positions_bed['seqname'].unique()
    
    strands = ("+", "-")
    
    all_matches = []
    
    for seqname in seqnames:
        for strand in strands:
            msg = "Finding matches for seqname: {}, strand: {}".format(seqname, strand)
            logger.debug(msg)
            
            m_positions_strand = positions_bed['strand'] == strand
            m_positions_seqname = positions_bed['seqname'] == seqname
            m_positions = m_positions_strand & m_positions_seqname
            
            m_intervals_strand = intervals_bed['strand'] == strand
            m_intervals_seqname = intervals_bed['seqname'] == seqname
            m_intervals = m_intervals_strand & m_intervals_seqname
            
            positions = positions_bed.loc[m_positions, 'start']
            interval_starts = intervals_bed.loc[m_intervals, 'start']
            interval_ends = intervals_bed.loc[m_intervals, 'end']
            
            interval_info = intervals_bed[m_intervals]
            position_info = positions_bed[m_positions]
            
            matches = get_position_intersections(positions, 
                                               interval_starts, 
                                               interval_ends, 
                                               interval_info=interval_info,
                                               position_info=position_info)
            
            all_matches.append(matches)
            
    all_matches = collection_utils.flatten_lists(all_matches)
    return all_matches

###
#   The following specialized functions find exact matches between two groups
#   of intervals.
#
#   TODO: it is unclear if this is more efficient than just passing 
#       min_overlap=1 to the generic overlap functions.
###

def get_exact_interval_matches(a_starts, a_ends, a_info, b_starts, b_ends, b_info):
    """ This function finds all of the intervals from a which exactly match
        with intervals from b. It allows arbitrary information to be attached
        with each of the intervals, and it returns a list of tuples of the
        information from each exact match.
        
        Internally, a variant of the "chrom-sweep" algorithm is used to
        efficiently find the intersections. See [Pedersen et al., Genome Biology
        17:118, 2016] for a good description of this algorithm.

        The complexity of the algorithm is O(m*k+n), where m and n are the
        number of intervals in a and b, and k is the maximum number of intervals
        which overlap.

        Args:

            a_start, a_ends, b_starts, b_ends (list-likes of ints): lists of the start 
                and end positions of the intervals. The x_start and x_end arrays should
                match, so, for example, a_starts[i] matches to a_ends[i].
        
            a_info, b_info (list-likes): a list of additional information associated
                with the intervals. The container must be index-able, and it should
                match with the interval_starts and interval_ends. So 
                a_info[i] should give the information associated with a_start[i].

            N.B. The intervals will be converted to numpy arrays and sorted as the 
            first step of this function, so they do not need to be pre-sorted.

        Returns:
            list of 2-tuples: a list of the additional information for each pair
                of intervals which exactly overlap.
    """
    # convert to numpy
    a_starts = np.array(a_starts)
    a_ends = np.array(a_ends)
    a_info = np.array(a_info)
    
    b_starts = np.array(b_starts)
    b_ends = np.array(b_ends)
    b_info = np.array(b_info)
    
    num_a_intervals = len(a_starts)
    num_b_intervals = len(b_starts)

    # sort by starts
    sorted_a_indices = np.argsort(a_starts)
    a_starts = a_starts[sorted_a_indices]
    a_ends = a_ends[sorted_a_indices]
    a_info = a_info[sorted_a_indices]
    
    sorted_b_indices = np.argsort(b_starts)
    b_starts = b_starts[sorted_b_indices]
    b_ends = b_ends[sorted_b_indices]
    b_info = b_info[sorted_b_indices]
    
    # bookkeeping information
    next_a_interval = 0
    next_a_start = a_starts[0]
    next_a_end = a_ends[0]

    next_b_interval = 0
    next_b_start = b_starts[0]
    next_b_end = b_ends[0]

    cache = []

    matches = []
    while next_a_interval < num_a_intervals:

        # get whichever interval comes next
        if next_a_start < next_b_start:

            # check if this exactly matches anything in the cache
            for c in cache:
                starts = b_starts[c] == next_a_start
                ends = b_ends[c] == next_a_end
                if starts and ends:
                    # we have a match
                    matches.append((a_info[next_a_interval], b_info[c]))

            # now, discard anything which starts before this
            cache = [c for c in cache if b_starts[c] >= next_a_start]
            
            # and advance
            next_a_interval += 1

            # if we checked the last a interval, then we will not
            # find any more matches
            if next_a_interval == num_a_intervals:
                break
            next_a_start = a_starts[next_a_interval]
            next_a_end = a_ends[next_a_interval]

        else:

            # just add it to the cache
            cache.append(next_b_interval)
            
            # and advance
            next_b_interval += 1

            if next_b_interval < num_b_intervals:
                next_b_start = b_starts[next_b_interval]
                next_b_end = b_ends[next_b_interval]
            else:
                next_b_start = np.inf
                
    return matches

def get_exact_block_matches(matches, block_counts_a, block_counts_b=None, 
        block_id_index=None):

    """ This function finds pairs of transcripts (or whatever outer-level
        object is considered) which have exact interval matches for all of
        their blocks (i.e., exons). Roughly, it does this by counting the
        number of exact matches for each pair of transcripts. Then, it
        checks if the number of matches is equal to the number of blocks
        in each transcript.

        For example, this will not consider transcripts with different
        number of blocks as ever being an exact match.

        The complexity of this algorithm is O(m), where m is the number
        of matches.

        Args:
            matches: a list of pairs of interval matches, presumably from 
                get_exact_interval_matches. Each element of the pair must
                contain some way to identify the transcript to which the
                transcript to which the interval belongs.

            block_counts_a (dict-like): the number of blocks (exons) in each
                transcript in a

            block_counts_b (dict-like or None): the number of blocks (exons) in
                each transcript in b, or None. If None, then block_counts_a is
                assumed to contain the counts for both a and b.

            block_id_index (int or None): the index of the transcript identifier
                within the match info objects, or None. If None, then the match
                info is assumed to be the identifier. These identifiers must
                match the keys in the block_counts dicts.

        Returns:
            list of 2-tuples: a list of the transcript identifers which are
                exact matches
    """
    # keep track of the exact matches between different transcripts
    match_count = collections.defaultdict(int)
    
    if block_counts_b is None:
        block_counts_b = block_counts_a

    for match in matches:
        a_match = match[0]
        b_match = match[1]
        
        if block_id_index is not None:
            a_match = a_match[block_id_index]
            b_match = b_match[block_id_index]

        key = (a_match, b_match)

        match_count[key] += 1
        
    # now, check if we see enough matches for any of the pairs
    exact_matches = []

    for (a_match, b_match), num_matches in match_count.items():
        num_transcript_exons = max(block_counts_a[a_match], block_counts_b[b_match])

        if num_matches == num_transcript_exons:
            match = (a_match, b_match)
            exact_matches.append(match)
            
    return exact_matches

def get_exact_bed_matches(bed_a, bed_b, seqname, strand):
    """ This function finds the exact transcript matches of two BED6+ data
        frames, on the indicated seqname and strand. This function is mostly
        a convenience wrapper around get_exact_block_matches.

        N.B. This function *does not* split BED12+ records. It asssumes they
            are already split, and the "id" field gives the transcript to which
            each block belongs.

        Args:
            bed_a (pd.DataFrame): a data frame corresponding to a BED6+ data set

            bed_b (pd.DataFrame): a data frame corresponding to a BED6+ data set

            seqname (str): the sequence on which matches will be looked

            strand (str): the strand on which sequences will be looked

        Returns:
            list of 2-tuples: a list of transcript pairs which are exact matches.
                The tuples give the "id" field of the matches
    """
    m_bed_a_seqname = bed_a['seqname'] == seqname
    m_bed_b_seqname = bed_b['seqname'] == seqname
    
    m_bed_a_strand = bed_a['strand'] == strand
    m_bed_b_strand = bed_b['strand'] == strand
    
    a_starts = bed_a.loc[m_bed_a_seqname & m_bed_a_strand, 'start']
    a_ends = bed_a.loc[m_bed_a_seqname & m_bed_a_strand, 'end']
    a_info = bed_a.loc[m_bed_a_seqname & m_bed_a_strand, 'id']
    
    b_starts = bed_b.loc[m_bed_b_seqname & m_bed_b_strand, 'start']
    b_ends = bed_b.loc[m_bed_b_seqname & m_bed_b_strand, 'end']
    b_info = bed_b.loc[m_bed_b_seqname & m_bed_b_strand, 'id']
    
    block_counts_a = bed_a[m_bed_a_seqname & m_bed_a_strand].groupby('id').size()
    block_counts_a = block_counts_a.to_dict()
    
    block_counts_b = bed_b[m_bed_b_seqname & m_bed_b_strand].groupby('id').size()
    block_counts_b = block_counts_b.to_dict()
    
    if len(a_starts) == 0 or len(b_starts) == 0:
        msg = ("seqname: {}, strand: {}. Did not find any transcripts for a or "
            "b. len(a_starts): {}, len(b_starts): {}".format(seqname, strand, 
            len(a_starts), len(b_starts)))
        logger.warning(msg)
        return []
    
    matches = get_exact_interval_matches(a_starts, a_ends, a_info, b_starts, b_ends, b_info)
    exact_block_matches = get_exact_block_matches(matches, block_counts_a, block_counts_b)
    
    return exact_block_matches

def get_all_exact_bed_matches(bed_a, bed_b):
    """ This function finds the exact transcript matches of two BED6+ data
        frames.

        N.B. This function *does not* split BED12+ records. It asssumes they
            are already split, and the "id" field gives the transcript to which
            each block belongs.

        Args:
            bed_a (pd.DataFrame): a data frame corresponding to a BED6+ data set

            bed_b (pd.DataFrame): a data frame corresponding to a BED6+ data set

        Returns:
            pd.DataFrame: a data frame with the following fields:

                match: a 2-tuple that gives a transcript pair which is an exact
                    match. The tuples give the "id" field of the matches.

                seqname: the seqname of that transcripts

                strand: the strand of the transcripts
    """
    seqnames = bed_a['seqname'].unique()
    strands = ("+", "-")
    
    all_matches = []
    
    for seqname in seqnames:
        for strand in strands:
            matches = get_exact_bed_matches(bed_a, bed_b, seqname, strand)
            
            df = pd.DataFrame()
            df['match'] = matches
            df['seqname'] = seqname
            df['strand'] = strand
            
            all_matches.append(df)
            
    all_matches = pd.concat(all_matches)
    
    return all_matches

###
#   The following helper function is used to pull out the appropriate exons
#   based on the bed identifiers, if the exons are given
###
def _get_bed_exons(bed_a, bed_b, exons, exons_a, exons_b):
    """ Pull the exons for the bed entries from the appropriate exons group.

    This function is not intended for external use.

    Args:
        bed_a (pd.DataFrame): a data frame corresponding to a BED6+ data set

        bed_b (pd.DataFrame): a data frame corresponding to a BED6+ data set
        
        exons (pd.DataFrame): a data frame corresponding to a BED6+ data
            set containing entries from both A and B. The 'id's of A and
            'B' must match to the 'id's of the exons.

        exons_{a,b} (pd.DataFrame): a data frame corresponding to a BED6+ data
            set containing exons from A (or B) only. The 'id's of A (or B) must 
            match to the 'id's of the exons.
            
    Returns:
        bed_{a,b}: either the original data frames, or the appropriate entries
            from the exon data frames.

    Raises:
        RuntimeError: if exactly one of exons_{a,b} is passed, or if both
            exons and at least one of exons_{a,b} is passed.

    """
    # check that the exons we have makes sense
    exons_a_passed = exons_a is not None
    exons_b_passed = exons_b is not None
    if exons is not None:
        if exons_a_passed or exons_b_passed:
            msg = ("[get_bed_overlaps]: only exons or exons_{a,b} may be "
                "provided, but not both.")
            raise RuntimeError(msg)

    only_a  = exons_a_passed and not exons_b_passed
    only_b = not exons_a_passed and exons_b_passed
    if only_a or only_b:
        msg = ("[get_bed_overlaps]: either both or neither of exons_{a,b} must "
            "be provided.")
        raise RuntimeError(msg)

    if exons is not None:
        ids_a = set(bed_a['id'])
        m_exons_a = exons['id'].isin(ids_a)
        bed_a = exons[m_exons_a]
        
        ids_b = set(bed_b['id'])
        m_exons_b = exons['id'].isin(ids_b)
        bed_b = exons[m_exons_b]

    elif exons_a_passed:
        ids_a = set(bed_a['id'])
        m_exons_a = exons_a['id'].isin(ids_a)
        bed_a = exons_a[m_exons_a]
        
        ids_b = set(bed_b['id'])
        m_exons_b = exons_b['id'].isin(ids_b)
        bed_b = exons_b[m_exons_b]

    return bed_a, bed_b



###
#   The following set of functions is used to find arbitrary intersections
#   among bed entries.
###
interval_overlap = collections.namedtuple(
    "interval_overlap",
    "a_info,b_info,overlap"
)

transcript_overlap = collections.namedtuple(
    "transcript_overlap",
    "a_info,b_info,overlap,a_fraction,b_fraction"
)

def get_interval_overlaps(a_starts, a_ends, a_info, b_starts, b_ends, b_info):
    """ This function finds all of the intervals from a which match (overlap)
        with intervals from b. It allows arbitrary information to be attached
        with each of the intervals, and it returns a list of tuples of the
        information from each match, as well as the length of the overlap.
        
        Internally, a variant of the "chrom-sweep" algorithm is used to
        efficiently find the intersections. See [Pedersen et al., Genome Biology
        17:118, 2016] for a good description of this algorithm.

        The complexity of the algorithm is O(m*k+n), where m and n are the
        number of intervals in a and b, and k is the maximum number of intervals
        which overlap.

        Args:

            a_start, a_ends, b_starts, b_ends (list-likes of ints): lists of the start 
                and end positions of the intervals. The x_start and x_end arrays should
                match, so, for example, a_starts[i] matches to a_ends[i].
        
            a_info, b_info (list-likes): a list of additional information associated
                with the intervals. The container must be index-able, and it should
                match with the interval_starts and interval_ends. So 
                a_info[i] should give the information associated with a_start[i].

            N.B. The intervals will be converted to numpy arrays and sorted as the 
            first step of this function, so they do not need to be pre-sorted.

        Returns:
            list of 3-namedtuples: a list of the additional information for each pair
                of intervals which overlap, as well as the length of the overlap. The
                fields of the namedtuple are:
                
                    a_info: the additional information from the "a" list
                    b_info: the additional information from the "b" list
                    overlap: the length of the overlap
    """
    
    # convert to numpy
    a_starts = np.array(a_starts)
    a_ends = np.array(a_ends)
    a_info = np.array(a_info)
    
    b_starts = np.array(b_starts)
    b_ends = np.array(b_ends)
    b_info = np.array(b_info)
    
    num_a_intervals = len(a_starts)
    num_b_intervals = len(b_starts)

    # sort by starts
    sorted_a_indices = np.argsort(a_starts)
    a_starts = a_starts[sorted_a_indices]
    a_ends = a_ends[sorted_a_indices]
    a_info = a_info[sorted_a_indices]
    
    sorted_b_indices = np.argsort(b_starts)
    b_starts = b_starts[sorted_b_indices]
    b_ends = b_ends[sorted_b_indices]
    b_info = b_info[sorted_b_indices]
    
    # bookkeeping information
    next_a_interval = 0
    next_a_start = a_starts[0]
    next_a_end = a_ends[0]

    next_b_interval = 0
    next_b_start = b_starts[0]
    next_b_end = b_ends[0]

    a_cache = []
    b_cache = []

    matches = []
    while (next_a_interval < num_a_intervals) or (len(a_cache) != 0):

        # get whichever interval comes next
        if next_a_start < next_b_start:
            
            # a is first
            
            # remove everything in the b_cache which ends before this starts
            # i.e., keep things which end after this begins
            b_cache = [b for b in b_cache if b_ends[b] > next_a_start]

            # make sure we do not have any overlaps longer than the interval
            a_len = next_a_end - next_a_start
            
            # so whatever is in b_cache overlaps this
            for b in b_cache:
                overlap = b_ends[b] - next_a_start

                # it could happen that "b" starts before "a" starts and ends
                # after "a" ends. So, make sure the overlap is not longer than "a"
                overlap = min(a_len, overlap)

                io = interval_overlap(a_info[next_a_interval], 
                                        b_info[b], 
                                        overlap)

                matches.append(io)
            
            # add this to the a_cache
            a_cache.append(next_a_interval)
            
            # and advance
            next_a_interval += 1

            # make sure we don't run off the edge
            if next_a_interval < num_a_intervals:
                next_a_start = a_starts[next_a_interval]
                next_a_end = a_ends[next_a_interval]
            else:
                next_a_start = np.inf

        else:
            # b is first
                
            # remove everything in the a_cache which ends before this begins
            # i.e., keep things which end after this begins
            a_cache = [a for a in a_cache if a_ends[a] > next_b_start]
            
            # make sure we do not have any overlaps longer than the interval
            b_len = next_b_end - next_b_start
            
            # whatever is in the a_cache overlaps this
            for a in a_cache:
                overlap = a_ends[a] - next_b_start
                
                # it could happen that "a" starts before "b" starts and ends
                # after "b" ends. So, make sure the overlap is not longer than "b"
                overlap = min(b_len, overlap)

                io = interval_overlap(a_info[a], 
                                        b_info[next_b_interval], 
                                        overlap)
        
                matches.append(io)
            
            # add this to the b_cache
            b_cache.append(next_b_interval)
            
            # and advance
            next_b_interval += 1

            if next_b_interval < num_b_intervals:
                next_b_start = b_starts[next_b_interval]
                next_b_end = b_ends[next_b_interval]
            else:
                next_b_start = np.inf

    return matches

def get_transcript_overlaps(interval_overlaps):

    """ This function finds pairs of transcripts (or whatever outer-level
        object is considered) which have interval matches across multiple
        blocks and counts the total overlap.

        The complexity of this algorithm is O(m), where m is the number
        of matches.

        Args:
            matches: a list of pairs of interval matches and overlaps, 
                presumably from get_interval_matches. The a_info and
                b_info elements of the tuple must identify the transcripts.

        Returns:
            dictionary: (a_info, b_info) -> overlap
            
                A map from the pairs of overlapping transcripts to the total
                overlap between them.
    """

    # keep track of the total overlap between each pair
    total_overlaps = collections.defaultdict(int)
    
    for overlap in interval_overlaps:
        a_match = overlap.a_info
        b_match = overlap.b_info
        
        key = (a_match, b_match)

        total_overlaps[key] += overlap.overlap
        
    return total_overlaps

def get_transcript_overlap_fractions(transcript_overlaps, a_lengths_map, b_lengths_map, 
                                     min_a_overlap=0, min_b_overlap=0):
    """ This function finds the fraction of each transcript covered by
        each overlap found. Essentially, it just divides the length of 
        each transcript by the length of the overlap.
        
        The complexity of this algorithm is O(o), where o is the number
        of overlaps.
        
        Args:
            transcript_overlaps: dict-like mapping of pairs of transcript 
                identifiers to overlap lengths. Presumably, this is the
                output of get_transcript_overlaps.
                
            a_lengths_map: dict-like mapping from the transcript identifers
                in a to their lengths
                
            b_lengths_map: dict-like mapping from the transcript identifiers
                in b to their lengths
                
            min_a_overlap (float in [0,1]): a minimum fraction required
                for overlap of the "a" transcript. A fraction of "0" is
                interpreted to mean "at least one bp overlap."
            
            min_b_overlap (float in [0,1]): a minimum fraction required
                for overlap of the "b" transcript. A fraction of "0" is
                interpreted to mean "at least one bp overlap."
                
        Returns:
            list of transcript_overlap tuples: a list of the overlapping
                pairs of transcripts, the length of the overlap, and the 
                fraction of each transcript covered by that overlap
    """
    # check the overlap fractions
    validation_utils.check_range(min_a_overlap, 0, 1, variable_name="min_a_overlap")
    validation_utils.check_range(min_b_overlap, 0, 1, variable_name="min_b_overlap")

    transcript_fractions = []
    for (a_info, b_info), overlap in transcript_overlaps.items():
        a_length = a_lengths_map[a_info]
        b_length = b_lengths_map[b_info]

        a_fraction = overlap / a_length
        b_fraction = overlap / b_length
        
        if a_fraction < min_a_overlap:
            continue
            
        if b_fraction < min_b_overlap:
            continue

        to = transcript_overlap(a_info,
                                b_info,
                                overlap,
                                a_fraction,
                                b_fraction
                               )

        transcript_fractions.append(to)
        
    return transcript_fractions

def get_bed_overlaps(bed_a, bed_b, min_a_overlap=0, min_b_overlap=0, exons=None,
            exons_a=None, exons_b=None):
    """ This function finds the transcript overlaps of two BED6+ data
        frames. This function is mostly a convenience wrapper around 
        get_interval_overlaps, get_transcript_overlaps, and 
        get_transcript_overlap_fractions.
        
        Optionally, a minimum overlap fraction can be specified for either
        or both sets of transcripts. Only overlaps which meet the specified
        criteria will be returned.

        N.B. This function *does not* split BED12+ records. It asssumes they
            are already split, and the "id" field gives the transcript to which
            each block belongs.
        
        *** Alternatively, if the exons are passed in, then the 'id's from A and
            B will be used to extract the relevant exons. So, this function will
            not perform the splitting, but if it is already done, this can be
            used as a convenience wrapper.

            As a second alternative, separate exons can be passed for A and B.
            Then ids from the respective bed data frames will be used to find
            the appropriate exons.

        Args:
            bed_a (pd.DataFrame): a data frame corresponding to a BED6+ data set

            bed_b (pd.DataFrame): a data frame corresponding to a BED6+ data set
            
            min_a_overlap (float in [0,1]): a minimum fraction required
                for overlap of the "a" transcript. A fraction of "0" is
                interpreted to mean "at least one bp overlap."
            
            min_b_overlap (float in [0,1]): a minimum fraction required
                for overlap of the "b" transcript. A fraction of "0" is
                interpreted to mean "at least one bp overlap."

            exons (pd.DataFrame): a data frame corresponding to a BED6+ data
                set containing entries from both A and B. The 'id's of A and
                'B' must match to the 'id's of the exons.

            exons_a (pd.DataFrame): a data frame corresponding to a BED6+ data
                set containing exons from A only. The 'id's of A must match to 
                the 'id's of the exons.
                
            exons_b (pd.DataFrame): a data frame corresponding to a BED6+ data
                set containing exons from B only. The 'id's of B must match to 
                the 'id's of the exons.

        Returns:
            list of transcript_overlap tuples: a list of the overlapping
                pairs of transcript "id"s, the length of the overlap, and the 
                fraction of each transcript covered by that overlap

        Raises:
            RuntimeError: if exactly one of exons_{a,b} is passed, or if both
                exons and at least one of exons_{a,b} is passed.
    """
    bed_a, bed_b = _get_bed_exons(bed_a, bed_b, exons, exons_a, exons_b)
    
    # check the overlap fractions
    validation_utils.check_range(min_a_overlap, 0, 1, variable_name="min_a_overlap")
    validation_utils.check_range(min_b_overlap, 0, 1, variable_name="min_b_overlap")
    
    # first, get the total lengths of all transcripts
    bed_a['length'] = bed_a['end'] -  bed_a['start']
    bed_b['length'] = bed_b['end'] -  bed_b['start']

    a_groups = bed_a.groupby('id')
    a_lengths = a_groups['length'].sum()
    a_lengths_map = a_lengths.to_dict()

    b_groups = bed_b.groupby('id')
    b_lengths = b_groups['length'].sum()
    b_lengths_map = b_lengths.to_dict()
    
    seqnames = bed_a['seqname'].unique()
    strands = ("+", "-")
    all_transcript_overlaps = []
    
    seqnames_str = ','.join(seqnames)
    msg = "Looking for overlaps in seqnames: {}".format(seqnames_str)
    logger.debug(msg)
    
    for seqname in seqnames:
        for strand in strands:

            m_bed_a_seqname = bed_a['seqname'] == seqname
            m_bed_b_seqname = bed_b['seqname'] == seqname

            m_bed_a_strand = bed_a['strand'] == strand
            m_bed_b_strand = bed_b['strand'] == strand
            
            m_bed_a = m_bed_a_seqname & m_bed_a_strand
            m_bed_b = m_bed_b_seqname & m_bed_b_strand

            a_starts = bed_a.loc[m_bed_a, 'start']
            a_ends = bed_a.loc[m_bed_a, 'end']
            a_info = bed_a.loc[m_bed_a, 'id']

            b_starts = bed_b.loc[m_bed_b, 'start']
            b_ends = bed_b.loc[m_bed_b, 'end']
            b_info = bed_b.loc[m_bed_b, 'id']

            if len(a_starts) == 0 or len(b_starts) == 0:
                msg = ("seqname: {}, strand: {}. Did not find any transcripts for a or "
                    "b. len(a_starts): {}, len(b_starts): {}".format(seqname, strand, 
                    len(a_starts), len(b_starts)))
                logger.warning(msg)
                continue

            overlaps = get_interval_overlaps(a_starts, a_ends, a_info, b_starts, b_ends, b_info)
            transcript_overlaps = get_transcript_overlaps(overlaps)
            
            transcript_fractions = get_transcript_overlap_fractions(transcript_overlaps, 
                                                                    a_lengths_map, 
                                                                    b_lengths_map,
                                                                    min_a_overlap=min_a_overlap,
                                                                    min_b_overlap=min_b_overlap)

            all_transcript_overlaps.extend(transcript_fractions)
    
    return all_transcript_overlaps

###
#   These function handle various sorts of overlap operations, such as 
#   subtraction. Typically, they are just wrappers around get_bed_overlaps
#   with some pre- or post-processing.
###
def subtract_bed(bed_a, bed_b, min_a_overlap=0, min_b_overlap=0, exons=None,
            exons_a=None, exons_b=None):
    """ This function removes all entries from A which overlap any of the
        entries in B. Optionally, some minimum overlap can be given, in
        terms of fraction of overlap.

        N.B. This function *does not* split BED12+ records. It asssumes they
            are already split, and the "id" field gives the transcript to which
            each block belongs.

        *** Alternatively, if the exons are passed in, then the 'id's from A and
            B will be used to extract the relevant exons. So, this function will
            not perform the splitting, but if it is already done, this can be
            used as a convenience wrapper.

        Args:
            bed_a (pd.DataFrame): a data frame corresponding to a BED6+ data set

            bed_b (pd.DataFrame): a data frame corresponding to a BED6+ data set
            
            min_a_overlap (float in [0,1]): a minimum fraction required
                for overlap of the "a" transcript. A fraction of "0" is
                interpreted to mean "at least one bp overlap."
            
            min_b_overlap (float in [0,1]): a minimum fraction required
                for overlap of the "b" transcript. A fraction of "0" is
                interpreted to mean "at least one bp overlap."

            exons (pd.DataFrame): a data frame corresponding to a BED6+ data
                set containing entries from both A and B. The 'id's of A and
                'B' must match to the 'id's of the exons.
            
            exons_a (pd.DataFrame): a data frame corresponding to a BED6+ data
                set containing exons from A only. The 'id's of A must match to 
                the 'id's of the exons.
                
            exons_b (pd.DataFrame): a data frame corresponding to a BED6+ data
                set containing exons from B only. The 'id's of B must match to 
                the 'id's of the exons.


        Returns:
            set of transcript ids: a set of transcript identifiers (from the 
                "id" column) of the A set which do not overlap the B set,
                respecting the specified overlaps.
    """
        
    # first, find the overlaps
    overlaps = get_bed_overlaps(bed_a, bed_b, 
                                min_a_overlap=min_a_overlap, 
                                min_b_overlap=min_b_overlap,
                                exons=exons, exons_a=exons_a, exons_b=exons_b)
    
    # now, remove all of the A's which had an overlap
    a_no_match = set(bed_a['id'])
    for overlap in overlaps:
        a_id = overlap.a_info
        a_no_match.discard(a_id)
        
    # and return everything that is left
    return a_no_match

###
#   The following functions are used to extract sequences in A which have
#   overlaps in B that are not necessarily exactly in the body of A
###
def get_entries_with_upstream_overlaps(
        bed_a, 
        bed_b, 
        upstream_window, 
        allow_overlaps=False, 
        exons=None,
        exons_a=None,
        exons_b=None
    ):

    """ This function finds all intervals of A which have upstream intervals of B.
        It always takes the strand of the intervals into account. By default, the
        function looks for intervals in B which are strictly upstream of the 
        intervals in B; optionally, overlaps of the A and B features are allowed.

        N.B. This function *does not* split BED12+ records. It asssumes they
            are already split, and the "id" field gives the transcript to which
            each block belongs.

        *** Alternatively, if the exons are passed in, then the 'id's from A and
            B will be used to extract the relevant exons. So, this function will
            not perform the splitting, but if it is already done, this can be
            used as a convenience wrapper.

        Args:
            bed_a (pd.DataFrame): a data frame corresponding to a BED6+ data set

            bed_b (pd.DataFrame): a data frame corresponding to a BED6+ data set

            upstream_window (int): how far upstream of intervals in A to look for
                intervals in B.
            
            allow_overlaps (bool): whether to require intervals in B are strictly
                upstream of the intervals in A, or are allowed to overlap.

            exons (pd.DataFrame): a data frame corresponding to a BED6+ data
                set containing entries from both A and B. The 'id's of A and
                'B' must match to the 'id's of the exons.
            
            exons_a (pd.DataFrame): a data frame corresponding to a BED6+ data
                set containing exons from A only. The 'id's of A must match to 
                the 'id's of the exons.
                
            exons_b (pd.DataFrame): a data frame corresponding to a BED6+ data
                set containing exons from B only. The 'id's of B must match to 
                the 'id's of the exons.

        Returns:
            list of transcript_overlap tuples: a list of the overlapping
                pairs of transcript "id"s, the length of the overlap, and the 
                fraction of each transcript covered by that overlap

            N.B. The length and fraction of the overlap are not generally very
                meaningful here.
    """
    bed_a, bed_b = _get_bed_exons(bed_a, bed_b, exons, exons_a, exons_b)
    
    # we will always need strand information about A    
    m_a_forward = bed_a['strand'] == '+'
    m_a_reverse = bed_a['strand'] == '-'
        
    # check if we are given exons
    given_a_exons = (exons is not None) or (exons_a is not None)
    
    # if so, then we need to pull out only the first for each id in bed_a
    if given_a_exons:
        sorted_bed_a = bed_a.sort_values("start")

        # for forward-strand entries, the exon with the smallest start is first
        m_a_forward = sorted_bed_a['strand'] == '+'

        forward_groups = sorted_bed_a[m_a_forward].groupby("id", as_index=False)
        bed_a_upstream_forward = forward_groups.first()
        bed_a_upstream_forward = bed_a_upstream_forward.copy()

        # and vice verse for reverse-strand entries
        m_a_reverse = sorted_bed_a['strand'] == '-'

        reverse_groups = sorted_bed_a[m_a_reverse].groupby("id", as_index=False)
        bed_a_upstream_reverse = reverse_groups.last()
        bed_a_upstream_reverse = bed_a_upstream_reverse.copy()
               
    else:
        # otherwise, just use all of the A entries
        bed_a_upstream_forward = bed_a[m_a_forward].copy()
        bed_a_upstream_reverse = bed_a[m_a_reverse].copy()
    
    
    # first, we need to create regions of A which includes the relevant "upstream" regions
    bed_a_upstream_forward['end'] = bed_a_upstream_forward['start']
    bed_a_upstream_forward['start'] = bed_a_upstream_forward['start'] - upstream_window

    bed_a_upstream_reverse['start'] = bed_a_upstream_reverse['end']
    bed_a_upstream_reverse['end'] = bed_a_upstream_reverse['end'] + upstream_window

    bed_a_upstream = pd.concat([bed_a_upstream_forward, bed_a_upstream_reverse])
    
    upstream_overlaps = get_bed_overlaps(bed_a_upstream, bed_b)
    upstream_overlaps = {
        (o.a_info, o.b_info): o for o in upstream_overlaps
    }
    
    # check if we should filter the matches which extend into the "real" region
    a_ids = []
    b_ids = []

    if not allow_overlaps:
        a_ids = {o.a_info for o in upstream_overlaps.values()}
        b_ids = {o.b_info for o in upstream_overlaps.values()}
        
        # now, figure out whether to use exons
        if exons is None:
            m_a_id = bed_a['id'].isin(a_ids)
            m_b_id = bed_b['id'].isin(b_ids)

            bed_a_filtered = bed_a[m_a_id]
            bed_b_filtered = bed_b[m_b_id]
            
        else:            
            m_a_id = exons['id'].isin(a_ids)
            m_b_id = exons['id'].isin(b_ids)

            bed_a_filtered = exons[m_a_id]
            bed_b_filtered = exons[m_b_id]
            
        overlaps = get_bed_overlaps(bed_a_filtered, bed_b_filtered)
        
        # remove overlapping tuples
        for o in overlaps:
            upstream_overlaps.pop((o.a_info, o.b_info), None)
        
    res = list(upstream_overlaps.values())
    res = sorted(res)

    return res


###
#   The following set of functions is used to extract spliced sequences of
#   bed entries from fasta files.
###

def get_bed_sequence(bed_entry, seq_sequence, split_exons=True):
    """ This function finds the (spliced) sequence for the given bed entry
        from the given seq_sequence. It *does not* (and cannot) check that 
        seq_sequence matches the seqname of the provided bed entry.

        Args:
            bed_entry (dict-like): a dictionary-like object (such as a pandas
                Series) which uses the bed12_field_names as fields.

            seq_sequence (string): a string which contains the "chromosome"
                sequence to which the bed_entry belongs. This function *does
                not* verify that this is the correct sequence for the bed
                entry.

            split_exons (bool): whether to split the exons of the bed entry

        Returns:
            2-tuple:
                string: the header, which is just the "id" field of the bed entry

                string: the sequence for this entry, including splicing if 
                    specified
    """
    genomic_start = bed_entry['start']
    header = bed_entry['id']
    strand = bed_entry['strand']

    if not split_exons:
        genomic_end = bed_entry['end']
        transcript_sequence = seq_sequence[genomic_start:genomic_end]

    else:
        
        exon_starts = np.fromstring(bed_entry['exon_genomic_relative_starts'], 
                                        sep=',', dtype=int)

        exon_lengths = np.fromstring(bed_entry['exon_lengths'], 
                                        sep=',', dtype=int)
        
        genomic_starts = exon_starts + genomic_start
        genomic_ends = genomic_starts + exon_lengths

        transcript_sequence = [ seq_sequence[s:e] 
                                    for s,e in zip(genomic_starts, genomic_ends)]

        transcript_sequence = ''.join(transcript_sequence)

    if strand == '-':
        transcript_sequence = reverse_complement(transcript_sequence)
        
    return (header, transcript_sequence)

def get_all_bed_sequences(bed, fasta_file, bed_4=False, split_exons=True, progress_bar=True):
    """ This function extracts all of the sequences for entries in the given
        bed file. It is mostly a wrapper around get_bed_sequence, so please see
        that function for more information.

        Args:
            bed (string or pd.DataFrame): either a path to a bed file or a data
                frame containing bed entries

            fasta_file (string): the path to the fasta file containing the
                chromosome sequences

            bed_4 : bool
                Whether to treat the file as a BED4 or BED12 file

            split_exons (bool): whether to split the exons of the bed entries

        Returns:
            list of 2-tuples:
                string: the header, which is just the "id" field of the bed entry

                string: the sequence for this entry, including splicing if 
                    specified

        raises:
            ValueError: if the number of columns does not work
    """
    if isinstance(bed, str):
        msg = "Reading bed file"
        logger.debug(msg)

        bed = read_bed(bed)

    # make sure we have enough columns
    num_columns = len(bed.columns)
    if num_columns < 4:
        msg = "The bed file must contain at least four columns."
        raise ValueError(msg)
        
    if (num_columns < 12) and (not bed_4):
        msg = ("The bed file must contain at least twelve columns unless the "
            "--bed4 flag is given.")
        raise ValueError(msg)
    
    msg = "Opening fasta file"
    logger.debug(msg)

    fasta = fastx_utils.get_read_iterator(fasta_file)

    msg = "Finding bed entry sequences"
    logger.debug(msg)

    all_transcript_sequences = []


    for (seqname, sequence) in fasta:
        msg = "Processing seqname: {}".format(seqname)
        logger.debug(msg)
        
        seqname = seqname.split(" ")[0]
        m_seqname = bed['seqname'] == seqname

        if len(bed[m_seqname]) == 0:
            continue
        
        transcript_sequences = parallel.apply_df_simple(bed[m_seqname], 
                                                    get_bed_sequence, sequence, split_exons, 
                                                    progress_bar=progress_bar)
        
        all_transcript_sequences.extend(transcript_sequences)

    msg = "Found all bed entry sequences"
    logger.debug(msg)

    return all_transcript_sequences

###
#   The following functions are for high-level operations on bed files.
###
    
def sort(bed, seqname_order=None, transcript_ids=None):
    """ This file sorts the given BED data frame. If sorts the entries by seqname,
        start and end. Optionally, an order can be given which gives the precedence
        by which the seqnames should be sorted. Otherwise, they are sorted
        alphabetically.

        In the case of ties among seqname, start and end, the entries are sorted
        according to their order in the given data frame. Optionally, a list of
        'id's can be given, and the ties will be resolved according to the 
        precedence in the list.

        This order for sorting is chosen because it is the same order used by STAR.

        Args:
            bed (string or pd.DataFrame): either a path to a bed file or a data
                frame containing bed entries. If bed is a data frame, it must
                follow the bed12_field_name naming conventions.

            seqname_order (string or np.array): If this file is given, then the 
                bed entries will be sorted according to the order of seqnames in 
                this file. The file should contain one seqname on each line in
                the file. The file can also contain lines which do not include a
                seqname. They will not affect sorting.

            transcript_ids (np.array): If this array is given, then ties will
                broken according to the order of the transcripts in this list

        Returns:
            pd.DataFrame: a sorted data frame
    """
    
    # first, check if we need to read a file
    if isinstance(bed, str):
        msg = "Reading bed file"
        logger.debug(msg)

        bed12_df = read_bed(bed)
    else:
        bed12_df = bed.copy()


    # extract the transcript ids for breaking ties
    if transcript_ids is None:
        transcript_ids = np.array(bed12_df['id'])

    _, idx = np.unique(transcript_ids, return_index=True)
    unique_transcript_ids = transcript_ids[np.sort(idx)]
    id_map = {tid: i for i, tid in enumerate(unique_transcript_ids)}
    bed12_df['id_index'] = bed12_df['id'].map(id_map)

    # if the chr_name_file is given, also use it for sorting
    if seqname_order is not None:
        if isinstance(seqname_order, str):
            if os.path.exists(seqname_order):
                msg = "Using chr_name_file: {}".format(seqname_order)
                logger.debug(msg)

                seqname_order = pd.read_csv(seqname_order, header=None)
                seqname_order = np.array(seqname_order[0])

            else:
                msg = ("The chr_name_file was specified, but the file does not "
                    "exist: {}".format(seqname_order))
                raise FileNotFoundError(msg)

        seqname_order_map = {c: i for i, c in enumerate(seqname_order)}
        bed12_df['seqname_index'] = bed12_df['seqname'].map(seqname_order_map)
    else:
        bed12_df['seqname_index'] = bed12_df['seqname']

    # now do the actual sorting
    sort_fields = ['seqname_index', 'start', 'end', 'id_index']
    bed12_df = bed12_df.sort_values(sort_fields)

    return bed12_df

def concatenate(bed_list, sort_bed=False, out=None):
    """ Concatenate a list of bed files or data frames into a single bed-style
    data frame. Optionally, the concatenated data frame can be sorted and/or
    written to disk.

    Parameters
    ----------
    bed_list: list of strings and/or bed-like data frames
        The bed files or data frames to concatenate. All of the string entries
        will be treated as filenames, and the respective file contents will be
        included in the concatenated list.

    sort_bed: bool
        Whether to sort the concatenated data frame

    out: string or None
        If a string is given, the concatenated bed file will be written to the
        specified file. If specified, the bed records are sorted before being 
        written to disk. The (sorted) concatenated data frame is still returned.

    Returns
    -------
    concatenated_bed: bed-style data frame
        The concatenation of all of the bed records in the list
    """
    bed_list = [get_bed_df(b) for b in bed_list]
    concatenated_bed = pd.concat(bed_list)

    if sort_bed:
        concatenated_bed = sort(concatenated_bed)

    if out is not None:
        write_bed(concatenated_bed, out)

    return concatenated_bed

###
#   Utility functions for converting other formats to bed
###
def read_bam_as_bed(bam, progress_bar=True, logger=logger):
    """ Reads alignments from a bam file into a bed6+1 data frame. The
    additional field is the length of the read.

    Parameters
    ----------
    bam: string or pysam.AlignmentFile
        Either the path to a bam file or an open pysam.AlignmentFile. See
        bam_utils.get_pysam_alignment_file for more details.

    progress_bar: bool
        Whether to show a tqdm progress bar
    
    logger: logging.Logger
        The logger to use

    Returns
    -------
    bam_df: bed6+1 style pd.DataFrame
        The alignments as a bed data frame. The extra field is the read length.
    """
    # first, make sure we have an alignment file
    bam = bam_utils.get_pysam_alignment_file(bam)

    msg = "Counting the number of alignments"
    logger.debug(msg)
    num_alignments = bam.count()

    alignments = bam.fetch()

    lengths = np.zeros(num_alignments, dtype=int)
    starts = np.zeros(num_alignments, dtype=int)
    ends = np.zeros(num_alignments, dtype=int)
    seqnames = np.full(num_alignments, None, dtype=object)
    strands = np.full(num_alignments, None, dtype=object)

    msg = "Extracting the alignment coordinates"
    logger.debug(msg)

    for i, a in enumerate(tqdm.tqdm(alignments, leave=True, file=sys.stdout, total=num_alignments)):
        starts[i] = a.reference_start
        ends[i] = a.reference_end

        lengths[i] = a.qlen
        seqnames[i] = a.reference_name
        strands[i] = "+"
        
        if a.is_reverse:
            strands[i] = "-"

    msg = "Constructing data frame from alignment coordinates"
    logger.debug(msg)

    alignment_df = pd.DataFrame()
    alignment_df['seqname'] = seqnames
    alignment_df['start'] = starts
    alignment_df['end'] = ends + 1
    alignment_df['id'] = range(len(alignment_df))
    alignment_df['score'] = 0
    alignment_df['strand'] = strands
    alignment_df['length'] = lengths

    return alignment_df