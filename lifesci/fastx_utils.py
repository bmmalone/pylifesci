""" Helpers for working with fasta and fastq sequence files
"""
import logging
logger = logging.getLogger(__name__)

import collections
import gzip
import io
import os
import re
import shutil
import sys
import tqdm

# biopython
import Bio.SeqIO.FastaIO
import Bio.SeqIO.QualityIO

import pyllars.pandas_utils as pd_utils
import pyllars.string_utils as string_utils
import pyllars.utils

def get_read_iterator(fastx_file, is_fasta=True):
    """ This function returns an iterator (technically, a generator) for the
        reads in the specified file. The function does not attempt to guess 
        the type of the file should be specified correctly.

        The function *does* guess if a file is gzipped or not. In particular,
        if the name of the file ends in "gz", it is treated as a gzipped file.
        Otherwise, it is treated as plain text.

        Finally, if the passed file is an instantiation of io.IOBase (such as
        an already-open file), then the function just uses the file handle.

        Args:
            fastx_file (string or file-like): the path to the fasta or fastq 
                file. Alternatively, this can be an open file handle.

            is_fasta (bool): whether the file is a fasta file or not. If this
                value is false, then the file is treated as a fastq file.

        Returns:
            generator: which (lazily) iterates over the reads in the file.
                In particular, for fasta files, a generator from the class
                Bio.SeqIO.FastaIO.SimpleFastaParser is returned, while for
                fastq files, a Bio.SeqIO.QualityIO.FastqGeneralIterator is
                returned. Please see the relevant documentation from biopython
                for exact details of these semantics.
    """

    if isinstance(fastx_file, io.IOBase):
        msg = "Gussing that fastx_file is an open file handle"
        logger.debug(msg)

        file_handle = fastx_file
        
    elif fastx_file.endswith("gz"):

        msg = "Guessing that fastx_file is a gzipped file"
        logger.debug(msg)

        file_handle = gzip.open(fastx_file, 'rt')

    else:
        msg = "Guessing that fastx_file is a plain text file"
        logger.debug(msg)

        file_handle = open(fastx_file, 'rU')

    if is_fasta:
        read_iterator = Bio.SeqIO.FastaIO.SimpleFastaParser(file_handle)
    else:
        read_iterator = Bio.SeqIO.QualityIO.FastqGeneralIterator(file_handle)

    return read_iterator

def get_length_distribution(fastx_file, is_fasta=True):
    """ Count the reads of each length in the given fasta file.

    Parameters
    ----------
    fastx_file: string or file handle
        A path or open file handle to the fastx file. In practice, this can
        be anything interpretable by get_read_iterator. Thus, if the path is
        to a gzipped file, nothing special needs to be done.

    is_fasta: bool
        True if the fastx_file is in fasta format (and False if it is fastq).

    Returns
    -------
    length_distribution_df: pd.DataFrame
        A data frame with the columns: length, count
    """
    fastx = get_read_iterator(fastx_file, is_fasta=is_fasta)

    length_distribution = collections.defaultdict(int)

    for read in fastx:
        read_len = len(read[1])
        length_distribution[read_len] += 1

    length_distribution_df = pd_utils.dict_to_dataframe(
        length_distribution,
        key_name='length',
        value_name='count'
    )

    return length_distribution_df

def get_read_count(filename, is_fasta=True):
    """ This function counts the number of reads in the specified file. The 
        function does not attempt to guess the type of the file should be 
        specified correctly.

        The function *does* guess if a file is gzipped or not. In particular,
        if the name of the file ends in "gz", it is treated as a gzipped file.
        Otherwise, it is treated as plain text.

        Finally, if the passed file is an instantiation of io.IOBase (such as
        an already-open file), then the function just uses the file handle.

        N.B. If a file handle is passed in, it will be "consumed" by this
            function, so the file handle would likely need to be reset.

        Args:
            filename (string): the path to the fasta or fastq file.

            is_fasta (bool): whether the file is a fasta file or not. If this
                value is false, then the file is treated as a fastq file.

        Returns:
            int: the number of reads in the file
    """
    read_iterator = get_read_iterator(filename, is_fasta)
    num_reads = sum(1 for read in read_iterator)
    return num_reads

def get_fasta_dict(filename, is_fasta=True, key_fn=None):
    if key_fn is None:
        key_fn = lambda x: x

    read_iterator = get_read_iterator(filename, is_fasta)
    read_dict = {
        key_fn(r[0]): r[1] for r in read_iterator
    }
    return read_dict


def get_fastq_qual_dict(filename, key_fn=None):
    if key_fn is None:
        key_fn = lambda x: x

    read_iterator = get_read_iterator(filename, is_fasta=False)
    qual_dict = {
        key_fn(r[0]): r[2] for r in read_iterator
    }
    return qual_dict


def remove_duplicate_sequences(fastas_in, fasta_out, compress=True, 
        lower_precedence_re=None, progress_bar=False):
    """ This function removes all duplicate sequences from a list of fasta 
        files and writes the remaining sequences back out as a fasta file.

        If desired, a regular expression can be given for "lower precedence"
        sequence identifiers. An example of using this precedence operator
        is in removing duplicate sequences from a fasta file which combines
        de novo assembled transcripts and annotated ones. In case a de novo
        transcript matches an annotated one, we would prefer to keep only
        the annotated transcript and identifier. Thus, we would pass an RE
        matching the de novo assembled identifiers (which have a lower
        precedence).

        If a precedence re is not given, or two identifiers have the same
        precedence, the first identifier encountered will be kept.

        This function is not implemented particularly efficiently.

        Args:
            fasta_in (list of string) : paths to fasta files which contains duplicates
            fasta_out (string) : path to the output file
            compress (bool) : whether to gzip the output
            lower_precendence_re (string) : a regular expression that matches
                the identifiers of lower precendence transcripts. (See the
                description for more details)
            progress_bar (bool) : whether to show a progress bar for reading
                the sequences and removing duplicates.

        Returns:
            None, but the output file is written
    """
    # we will find duplicates by storing the sequences as keys in a dictionary
    seq_dict = {}

    for fasta_in in fastas_in:

        seqs = get_read_iterator(fasta_in, is_fasta=True)

        if lower_precedence_re is not None:
            lp_re = re.compile(lower_precedence_re)
        else:
            lp_re = None

        iter_ = seqs

        if progress_bar:
            total = get_read_count(fasta_in, is_fasta=True)
            iter_ = tqdm.tqdm(seqs, leave=True, file=sys.stdout, total=total)
       
        for seq_id, seq in iter_:
            # if we have not seen this sequence before, just add it to the dictionary
            existing_id = seq_dict.get(seq, None)
            if existing_id is None:
                seq_dict[seq] = seq_id
            elif lp_re is not None:
                # check if one of the id's is of lower precedence

                # "match" returns None when there is no match
                new_is_low = (lp_re.match(seq_id) is None)
                existing_is_low = (lp_re.match(existing_id) is None)

                # we only care in the case that the new id is of higher precedence
                if existing_is_low and not new_is_low:
                    seq_dict[seq] = seq_id
            else:
                # otherwise, skip the duplicate
                continue

    id_dict = [(v,k) for k,v in seq_dict.items()]
    write_fasta(id_dict, fasta_out, compress=compress, wrap=True)

def _write_fasta_entry(out, header, seq, wrap=False):
    """ This helper function writes a single fasta entry to the given file.
    """
    out.write(">")
    out.write(header)
    out.write("\n")

    if wrap:
        seq = string_utils.simple_fill(seq)

    out.write(seq)
    out.write("\n")

def write_fasta(seqs, filename, compress=True, wrap=True, progress_bar=False):
    """ This function writes the provided sequences to a fasta file. The input
        is given as a list of tuples in which the first item is the fasta header 
        and the second is the sequence. The contents of the sequences are not 
        validated, so nucleotide, amino acid, nonstandard, etc., sequences are 
        all handled.

        Args:
            seqs (list of tuples) : a list of headers and sequences

            filename (string) : the name of the output file

            compress (boolean) : whether to gzip the output

            wrap (bool): whether to wrap the sequence to 80 characters

            progress_bar (boolean) : whether to show a progress bar

        Returns:
            None
    """
    if compress:
        out = gzip.open(filename, 'wt')
    else:
        out = open(filename, 'w')

    if progress_bar:
        seq_iter = tqdm.tqdm(seqs, leave=True, file=sys.stdout)
    else:
        seq_iter = seqs

    for header, seq in seq_iter:
        _write_fasta_entry(out, header, seq, wrap)
        
    out.close()

def _write_fastq_entry(out, header, seq, qual_scores, wrap=False):
    """ This helper function writes a single fastq entry to the given file.
    """
    out.write("@")
    out.write(header)
    out.write("\n")

    if wrap:
        seq = string_utils.simple_fill(seq)

    out.write(seq)
    out.write("\n")

    out.write("+")
    out.write(header)
    out.write("\n")

    if wrap:
        qual_scores = string_utils.simple_fill(qual_scores)

    out.write(qual_scores)
    out.write("\n")

def write_fastq(seqs, quals, filename, wrap=False, progress_bar=False):
    """ This function writes the provided sequences and qualities to a fastq file. The input
        is given as two dictionaries in which the keys in both are the fasta headers. The
        seqs dictionary should map from the header to the sequence, and the quals dictionary
        should map from the header to the quality scores. The contents of the sequences are 
        not validated, so nucleotide, amino acid, nonstandard, etc., sequences are all handled.
        
        N.B. All sequences in the seqs dictionary are written to the file; the quals dictionary
             must contain at least these headers. Other headers in the quals dictionary are
             just ignore (so the headers of quals can be a superset of seqs).

        Args:
            seqs (dictionary) : a mapping from the header to the sequences

            quals (dictionary) : a mapping from the header to the quality scores

            filename (string) : the name of the output file. compression will be guess based on the file extension

            wrap (bool): whether to wrap the sequence and quality scores to 80 chars

            progress_bar (boolean) : whether to show a progress bar

        Returns:
            None
    """

    out = pyllars.utils.open_file(filename, 'w')

    if progress_bar:
        seq_iter = tqdm.tqdm(seqs.items(), leave=True, file=sys.stdout)
    else:
        seq_iter = seqs.items()

    for header, seq in seq_iter:
        qual_scores = quals[header]
        _write_fastq_entry(out, header, seq, qual_scores)
        
    out.close()

def _check_fastq_read(read, read_index):
    """ This is a helper function for validing a single fastq read. It returns
        a list of all errors for this read. It is not intended for external use.
    """
    errors = []

    # the identifier line
    identifier_line_number = read_index * 4

    # the biopython iterator we use internally already takes care of
    # the sequence identifier
    #if len(read[0]) < 2:
    #    msg = ("ERROR on Line {}: The sequence identifier line was too "
    #        "short.".format(identifier_line_number))
    #    errors.append(msg)

    #if read[0][0] != "@":
    #    msg = ("ERROR on Line {}: First line of a sequence does not begin "
    #        "wtih @".format(identifier_line_number))
    #    errors.append(msg)

    #if read[0][1] != " ":
    #    msg = ("ERROR on Line {}: No Sequence Identifier specified before "
    #        "the comment.".format(identifier_line_number))
    #    errors.append(msg)

    # the raw sequence line
    sequence_line_number = read_index*4 + 1
    sequence_len = len(read[1])

    if sequence_len == 0:
        msg = ("ERROR on Line {}: Raw Sequence is shorter than the min read "
            "length: 0".format(sequence_line_number))
        errors.append(msg)

    # the plus line
    plus_line_number = read_index*4 + 2

    # the biopython iterator does not include the plus line
    # so skip this
    #if read[2][0] != '+':
    #    msg = ("ERROR on Line {}: Third line of a sequence does not begin "
    #        "with +".format(plus_line_number))
    #    errors.append(msg)

    # the quality string line
    quality_line_number = read_index*4 + 3

    # the quality string is at position 2 in the tuple (not 3 because the plus
    # line is omitted)
    quality_len = len(read[2])

    if quality_len != sequence_len:
        msg = ("ERROR on Line {}: Quality string length ({}) does not equal "
            "raw sequence length ({})".format(quality_line_number, quality_len, 
            sequence_len))
        errors.append(msg)

    return errors


def check_fastq_file(filename, break_on_error=True, raise_on_error=True, 
            logger=logger):

    """ This function checks that a fastq file is valid. Optionally, it 
        raises an exception if the file is invalid. Otherwise, it writes
        a "critical" warning message.

        The following criteria are used to determine if the fastq file is
        valid. They are a subset copied from: http://genome.sph.umich.edu/wiki/FastQ_Validation_Criteria

        Sequence Identifier Line:
            * Line is at least 2 characters long ('@' and at least 1 for the 
                sequence identifier)

            * Line starts with an '@'

            * Line does not contain a space between the '@' and the first 
                sequence identifier (which must be at least 1 character).

        Raw Sequence Line

            * A base sequence should have non-zero length.

        Plus Line

            * Must exist for every sequence.

            * Must start with a '+'

        Quality String Line

            * A quality string should be present for every base sequence.

            * Paired quality and base sequences should be of the same length.


        Args:
            filename (str): a path to the fastq file

            break_on_error (bool): whether to break on the first error found.
                This flag is only relevant if raise_on_error is False.

            raise_on_error (bool): whether to raise an OSError (if True) or log
                a "critical" message (if false)

            logger (logging.Logger): a logger for writing the message if an
                error is not raised

        Returns:
            bool: whether the file was valid

        Raises:
            OSError: if the file is not a valid fastq and raise_on_error is True
    """
    read_iterator = get_read_iterator(filename, is_fasta=False)

    is_valid = True

    # now, check the reads
    try:
        for i, read in enumerate(read_iterator):
            errors = _check_fastq_read(read, i)

            if len(errors) > 0:
                # we found some errors, so raise an error or break as specified
                is_valid = False

                msg = "The fastq file does not appear to be valid: {}".format(filename)
                errors.insert(0, msg)
                msg = '\n'.join(errors)

                if raise_on_error:
                    raise OSError(msg)

                logger.critical(msg)

                if break_on_error:
                    break

    except ValueError as e:
        is_valid = False

        msg = "The fastq file does not appear to be valid: {}".format(filename)
        msg = msg + "\n" + str(e)

        # the fastq iterator causes a ValueError if it runs into trouble
        if raise_on_error:
            raise OSError(msg)

        else:
            logger.critical(msg)

    # we have looked at either all the reads and found no errors, or found some
    # errors and quit as specified by the arguments
    return is_valid


def check_fasta_file(filename, break_on_error=True, raise_on_error=True, 
            logger=logger):
    """ This function checks the validity of a fasta file. This is somewhat
        difficult because the format is so loose; in contrast, for example,
        the length of the sequence and quality must match in fastq.

        This function just counts the sequences in the file.
        
        Args:
            filename (str): a path to the fasta file

            break_on_error (bool): whether to break on the first error found.
                This flag is only relevant if raise_on_error is False.

            raise_on_error (bool): whether to raise an OSError (if True) or log
                a "critical" message (if false)

            logger (logging.Logger): a logger for writing the message if an
                error is not raised

        Returns:
            bool: whether the file was valid

        Raises:
            OSError: if the file is not a valid fasta and raise_on_error is True
    """

    is_valid = True

    try:
        read_count = get_read_count(filename, is_fasta=True)

    except Exception as e:
        is_valid = False

        msg = "The fasta file does not appear to be valid: {}".format(filename)
        msg = msg + "\n" + str(e)

        if raise_on_error:
            raise OSError(msg)

        else:
            logger.critical(msg)

    # we have looked at either all the reads and found no errors, or found some
    # errors and quit as specified by the arguments
    return is_valid

