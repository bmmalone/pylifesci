"""
This module includes helpers for working with GTF files.

pyensembl includes many functions for querying, etc., GTF files. This
module does not replace that functionality. Rather, it adds light-weight
utilities for things like file conversion, etc.
"""
import logging
logger = logging.getLogger(__name__)

import csv
import numpy as np
import pandas as pd
import re
import shlex

import pyllars.pandas_utils as pd_utils
import lifesci.bed_utils as bed_utils

gtf_field_names = [
    "seqname",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attributes"
]

def read_gtf(filename, sep='\t', comment='#', field_names=None, 
    use_default_field_names=True):
    """ This function reads a GTF file into a pandas data frame. By default, it
        assumes comment lines start with a hash sign, and that the field names
        are omitted from the file.

        N.B. This function does not parse out attributes.

        This function forces the first field (seqname) to be treated as a string,
        but it does not affect the default types of the other fields.

        Args:
            filename (string) : the path to the GTF file

            sep, comments (string) : these arguments are passed through to
                pandas.read_csv

            use_default_field_names (bool) : if this is true, then the column
                names in lifesci.gtf_field_names will be used

            field_names (list of strings) : if use_default_field_names is false,
                then the columns in the data frame will have these names

        Returns:
            pandas.DataFrame: a data frame with the GTF records.
    """
    if use_default_field_names:
        gtf = pd.read_csv(filename, sep=sep, comment=comment, header=None,
            names=gtf_field_names)
    else:
        gtf = pd.read_csv(filename, sep=sep, comment=comment, header=None,
            names=field_names)

    seqname_field = gtf.columns[0]
    gtf[seqname_field] = gtf[seqname_field].astype(str)

    start_field = gtf.columns[3]
    end_field = gtf.columns[4]

    gtf[start_field] = gtf[start_field].astype(int)
    gtf[end_field] = gtf[end_field].astype(int)

    return gtf

def write_gtf(data_frame, filename, compress=True, use_default_fields=True, **kwargs):
    """ This function formats a data frame such that it will behave as
        a GTF file. In particular, it writes a tab-delimited file and prepends
        the hash mark (#) to the field names in the header so they are treated
        as comments by typical GTF parsers.

        The "start" and "end" features (4th and 5th columns) will be cast
        as integers.
    
        Args:
            data_frame (pandas.DataFrame) : a data frame representing bed objects

            filename (string) : the name of the output file

            compress (bool) : whether to gzip the output

            use_default_fields (bool) : whether to use all of the fields in
                the data frame or only those in gtf_field_names

            kwargs : these are passed through to the write_df function

        Returns:
            None
    """
    do_not_compress = not compress

    if use_default_fields:
        data_frame = data_frame[gtf_field_names]

    start_field = data_frame.columns[3]
    end_field = data_frame.columns[4]

    data_frame[start_field] = data_frame[start_field].astype(int)
    data_frame[end_field] = data_frame[end_field].astype(int)

    header = ['#{}'.format(c) for c in data_frame.columns]
    pd_utils.write_df(data_frame, filename, index=False, sep='\t', 
        header=header, do_not_compress=do_not_compress, quoting=csv.QUOTE_NONE, **kwargs)

def _get_gtf_entries(bed_entry, feature_type:str, source:str=None, 
        id_attribute:str="transcript_id"):
    """ Split the bed12 entry into multiple gtf entries, one for each block.
    
    The attributes for the gtf entries are taken as: "{field} {val}; ", where
    {col} is the name of the field and {val} is the respective value for all
    fields in bed_entry except for the standard bed12 fields (given by 
    lifesci.bed_utils.bed12_field_names).
    
    Additionally, if id_attribute is specified, then the "id" field in
    the bed entry will be used as the value for that attribute.
    
    Parameters
    ----------
    bed_entry: pd.Series
        The bed12 entry
        
    feature_type: string
        The "feature" field for the gtf entries (e.g., "exon" or "CDS")
        
    source: string or None
        The "source" field for the gtf entries. If None is given, "." is used
        for the source.
        
    id_attribute: string or None
        The name of the attribute for the "id" field. If None is given, then
        the "id" field from the bed entry is not used for anything.
        
    Returns
    -------
    gtf_entries: pd.DataFrame
        A data frame of the gtf entries for the bed entry blocks
    """
    # get the attributes

    # this comes from all of the "extra" columns, plus "id" for the 
    # transcript_id field
    attributes = bed_entry.drop(bed_utils.bed12_field_names)
    attributes = ["{} {}".format(k,shlex.quote(str(v))) 
                      for k,v in attributes.items()]
    
    if id_attribute is not None:
        attributes += ["{} {}".format(id_attribute, 
            shlex.quote(str(bed_entry['id'])))]

    attributes = "; ".join(attributes)
    
    if len(attributes) > 0:
        attributes += ";"

    # first, the exons
    start = bed_entry['start']
    block_lengths = bed_entry['exon_lengths'].split(",")
    block_rel_starts = bed_entry['exon_genomic_relative_starts'].split(",")

    block_lengths = np.array(block_lengths, dtype=int)
    block_rel_starts = np.array(block_rel_starts, dtype=int)

    block_starts = np.zeros(len(block_lengths), dtype=int)
    block_starts[1:] = np.cumsum(block_lengths)[:-1]

    exon_starts = np.array([
        bed_utils.get_gen_pos(rp, start, block_lengths, block_starts, block_rel_starts)
            for rp in block_starts
    ])

    exon_ends = exon_starts + block_lengths

    # move from half-open, base-0 (bed) to closed, base-1 (gtf)
    exon_starts += 1

    # there is no need to adjust exon stops, since the bed entry is open on 
    # that side
    
    if source is None:
        source = "."

    gtf_exons = [
        {
            "seqname": bed_entry['seqname'],
            "source": source,
            "feature": feature_type,
            "start": start,
            "end": end,
            "score": bed_entry['score'],
            "frame": ".",
            "strand": bed_entry['strand'],
            "attributes": attributes        
        } for start, end in zip(exon_starts, exon_ends)
    ]

    gtf_exons = pd.DataFrame(gtf_exons)
    return gtf_exons

def get_gtf_entries(bed_entry, source:str, id_attribute:str="transcript_id"):
    """ Split the bed12 entry into multiple gtf entries.
    
    "exon" entries will always be created for each block in the bed entry.
    Additionally, if the "thick_xxx" fields are given, then "CDS" entries
    will also be created.
    
    The attributes for the gtf entries are taken as: "{field} {val}; ", where
    {col} is the name of the field and {val} is the respective value for all
    fields in bed_entry except for the standard bed12 fields (given by 
    lifesci.bed_utils.bed12_field_names).
    
    Additionally, if id_attribute is specified, then the "id" field in
    the bed entry will be used as the value for that attribute.
    
    Parameters
    ----------
    bed_entry: pd.Series
        The bed12 entry
        
    feature_type: string
        The "feature" field for the gtf entries (e.g., "exon" or "CDS")
        
    source: string or None
        The "source" field for the gtf entries. If None is given, "."
        is used for the source.
        
    id_attribute: string or None
        The name of the attribute for the "id" field. If None is given,
        then the "id" field from the bed entry is not used for anything.
        
    Returns
    -------
    gtf_entries: pd.DataFrame
        A data frame of the gtf entries for the bed entry blocks. 
        
        * The entries are sorted according to "start".
        * The order of the columns is the same as gtf_utils.gtf_field_names
    """    
    gtf_exons = _get_gtf_entries(bed_entry, "exon", source)

    if bed_entry['thick_start'] > -1:
        bed_entry_cds = bed_utils.retain_thick_only(bed_entry)
        gtf_cds = _get_gtf_entries(bed_entry_cds, "CDS", source)

        gtf_exons = pd.concat([gtf_exons, gtf_cds])

    gtf_exons = gtf_exons.sort_values('start')
    gtf_exons = gtf_exons.reset_index(drop=True)
    return gtf_exons[gtf_field_names]

###
#   These functions are helpers for using GTF files
###
###
#
# Much of the parsing code is adapted from: https://gist.github.com/slowkow/8101481
#
###
import re
R_SEMICOLON = re.compile(r'\s*;\s*')
R_COMMA     = re.compile(r'\s*,\s*')
R_KEYVALUE  = re.compile(r'(\s+|\s*=\s*)')

def _get_gtf_value(value):
    if not value:
        return None

    # Strip double and single quotes.
    value = value.strip('"\'')

    # Return a list if the value has a comma.
    if ',' in value:
        value = re.split(R_COMMA, value)
    # These values are equivalent to None.
    elif value in ['', '.', 'NA']:
        return None

    return value

def parse_gtf_attributes(row):
    """ This function parses the attributes of a GTF entry, where the GTF fields
        are consistent with gtf_field_names. This returns a new GTF entry with
        the attributes added as keys.

        For "flag" attributes, the key and the value are both the name of the
        flag.

        Args:
            row (dict-like): a dictionary-like data structure (such as a dict 
                or pd.Series) which includes the gtf_field_names. Additionaly,
                the object must have a copy() method.

        Returns:
            :
            dict-like: an object of the same types as the input which includes
                all of the fields present in the input, as well as each
                attribute as a new field (e.g., key in a dictionary)
    """
    attributes = row['attributes']
    
    attributes = [x for x in re.split(R_SEMICOLON, attributes) if x.strip()]
    result = row.copy()

    for i, attribute in enumerate(attributes, 1):
        # It should be key="value".
        try:
            key, _, value = re.split(R_KEYVALUE, attribute, 1)
        # But sometimes it is just "value".
        except ValueError:
            key = 'INFO{}'.format(i)
            value = attribute
        # Ignore the field if there is no value.
        if value:
            result[key] = _get_gtf_value(value)
   
    return result


def _parse_gtf_group(rows):
    """ This is a helper function for parsing GTF attributes from a data frame.
        It is not intended for external use.
    """
    res = parallel.apply_df_simple(rows, parse_gtf_attributes)
    res = pd.DataFrame(res)
    return res


def parse_all_gtf_attributes(gtf_df, num_cpus=1, progress_bar=False, num_groups=100):
    """ This function parses the attributes of each entry in the given GTF
        data frame. It returns a new data frame with the attributes added
        as fields.

        See parse_gtf_attributes for more details.

        Args:
            gtf_df (pd.DataFrame): a data frame containing fields consistent
                with gtf_field_names

            num_cpus (int): the number of CPUs to use for extracting the features

            progress_bar (bool): whether to show a progress bar

            num_groups (int): the number of groups into which the data frame
                will be split for parallelization

        Returns:
            pd.DataFrame: a new data frame with each attribute added as a column
    """ 
    gtf = parallel.apply_parallel_split(
        gtf_df, 
        num_cpus, 
        _parse_gtf_group, 
        progress_bar=progress_bar, 
        num_groups=num_groups)

    gtf = pd.concat(gtf)
    return gtf