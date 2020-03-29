"""
Helper functions for working with the MEME suite: http://meme-suite.org/
"""
import logging
logger = logging.getLogger(__name__)

import collections
import numpy as np
import pandas as pd

import pyllars.collection_utils as collection_utils

meme_motifs = collections.namedtuple(
    "meme_motifs",
    "version,alphabet,strands,background_frequencies,motifs"
)

motif = collections.namedtuple(
    "motif",
    ("identifier,alternate_name,probability_matrix_attributes,"
     "probability_matrix,url")
)

def parse_meme_motifs(meme_motifs_file):
    """ Parse a meme motifs file into a list of motifs (and some metadata).
    
    The file format is described here: 
        http://meme-suite.org/doc/meme-format.html
    
    N.B. This implementation is likely quite fragile; it has only been tested
    with a few of the files distributed in the official release, version 4.
    
    In general, the field names, etc., follow directly from the official
    documentation, so please refer to it for more details.
    
    Parameters
    ----------
    meme_motifs_file: string
        The path to the meme motifs file
        
    Returns
    -------
    meme_motifs: a meme_motifs namedtuple.
        The fields of the named tuple are:
        - version. float. the version of the meme file
        
        - alphabet. string. the alphabet, as a single string
        
        - strands. list of strings. the specified strands,
            or both strands if none are specified
        
        - background_frequencies. dict. a map from each
            alphabet character to its frequency, or None if
            the background frequencies were not specified.
        
        - motifs. dict of motif identifiers to motif 
            namedtuples. the tuples have these fields.
            
            - identifier
            - alternate_name
            - probability_matrix_attributes
            - probability_matrix
            - url
    """
    version = None
    alphabet = None
    strands = ["+", "-"]
    background_frequencies = None
    motifs = {}
    
    read_background = False

    with open(meme_motifs_file) as f:

        for line in f:
            ll = line.strip().lower()

            # check for the version
            if ll.startswith("meme version"):
                s = ll.split()
                version = float(s[2])

            elif ll.startswith("alphabet"):
                alphabet = line.strip().split()[1]

            elif ll.startswith("end alphabet"):
                msg = "This parser does not support custom alphabet definition"
                raise SyntaxWarning(msg)

            elif ll.startswith("strands"):
                strands = ll.split()[1:]

            elif ll.startswith("motif"):
                # once we see motifs, break out
                next_motif_line = line.strip()
                break

            elif ll.startswith("background letter frequencies"):
                read_background = True
                background_frequencies = {}

            elif read_background:
                if len(ll) == 0:
                    read_background = False
                else:
                    s = line.strip().split()
                    b = collection_utils.list_to_dict(s, float)
                    background_frequencies.update(b)

        # keep parsing until we find no more motifs
        while len(next_motif_line) > 0:
            m, next_motif_line = _parse_meme_motif(next_motif_line, f)
            motifs[m.identifier] = m

    mm = meme_motifs(
        version=version,
        alphabet=alphabet,
        strands=strands,
        background_frequencies=background_frequencies,
        motifs=motifs
    )

    return mm

def _parse_meme_motif(motif_line, f):
    """ Parse the next meme motif from the file.

    Parameters
    ----------
    motif_line: string
        the motif name line

    f: file
        the open file, which should be pointing to the line immediately after
        the motif name line

    Returns
    -------
    motif: motif named tuple
        the motif structure

    next_motif_line, or "": string
        the next motif name line, or the empty string if this was the last
        motif in the file
    """
    s = motif_line.split()
    
    if len(s) == 2:
        motif_kw, identifier = s
    if len(s) == 3:
        motif_kw, identifier, alternate_name = s
    
    probability_matrix_attributes = None
    read_letter_probability = False
    all_letter_probabilities = []
    url = None
    
    for line in f:
        ll = line.strip().lower()

        if ll.startswith("letter-probability"):
            read_letter_probability = True

            s = line.split()
            pma = collection_utils.list_to_dict(s[2:])
            probability_matrix_attributes = pma

        elif ll.startswith("url"):
            read_letter_probability = False
            url = line.split()[1]

        elif read_letter_probability and len(ll) > 0:
            s = ll.split()
            letter_probabilities = [float(f) for f in s]
            all_letter_probabilities.append(letter_probabilities)
            
        elif ll.startswith("motif"):
            break
            
    letter_probability_matrix = np.array(all_letter_probabilities)
    
    m = motif(
        identifier=identifier,
        alternate_name=alternate_name,
        probability_matrix_attributes=probability_matrix_attributes,
        probability_matrix=letter_probability_matrix,
        url=url
    )
    
    return m, line.strip()
    

def read_fimo(fimo_file):
    """ This function reads the output of fimo into a data frame.

    Please see the documentation for more about fimo: 
        http://meme-suite.org/doc/fimo.html

    Parameters
    ----------
    fimo_file : string
        The path to the output of fimo

    Returns
    -------
    fimo_df : pd.DataFrame
        A data frame containing the following columns:
            pattern name : the name of the motif
            sequence name : the name of the sequence containing the motif
            start : the position where the motif starts
            stop : the position where the motif ends
            strand : the strand of the sequence
            score : see documentation
            p-value : the significance of the match
            q-value : the lowest FDR at which the match is significant
            matched sequence : the sequence which matched the motif
    """
    df = pd.read_csv(fimo_file, sep='\t')
    df.columns = [c.replace("#", "") for c in df.columns]
    return df

def _parse_ame_ranksum_line(line):
    if line.find("Ranksum p-values") == -1:
        return None
        
    # manually parse out the lines... ugg...
    s = line.strip().split()
    res = {
        "motif_name": s[5],
        "num_seqs": int(s[7]),
        "left_pvalue": float(s[10]),
        "right_pvalue": float(s[11]),
        "twotailed_pvalue": float(s[12]),
        "u_value": float(s[14]),
        "corrected_left_pvalue": float(s[17]),
        "corrected_right_pvalue": float(s[18]),
        "corrected_twotailed_pvalue": s[19]
    }

    # remove the final ")" from the corrected_twotailed_pvalue
    res["corrected_twotailed_pvalue"] = res["corrected_twotailed_pvalue"][:-1]
    res["corrected_twotailed_pvalue"] = float(res["corrected_twotailed_pvalue"])

    return res

def _parse_ame_fisher_line(line):
    # example line (split on spaces): 
    #   
    # [0] 1.
    # [1] Fisher's
    # [2] exact
    # [3] test
    # [4] p-value
    # [5] of
    # [6] motif
    # [7] M151_0.6
    # [8] (Hnrnph2)_(Homo_sapiens)_(RBD_0.99)
    # [9] top
    # [10] 141
    # [11] seqs:
    # [12] 5.872e-05
    # [13] (Corrected
    # [14] p-value:
    # [15] 0.01162)

    if line.find("Fisher's") == -1:
        return None

    s = line.strip().split()

    pvalue = float(s[12])
    corrected_pvalue = s[15]
    corrected_pvalue = corrected_pvalue[:-1]
    corrected_pvalue = float(corrected_pvalue)

    res = {
        "motif_name": s[7],
        "num_seqs": int(s[10]),
        "left_pvalue": pvalue,
        "right_pvalue": pvalue,
        "twotailed_pvalue": pvalue,
        "u_value": -1,
        "corrected_left_pvalue": corrected_pvalue,
        "corrected_right_pvalue": corrected_pvalue,
        "corrected_twotailed_pvalue": corrected_pvalue
    }
    return res


AME_METHODS = [
    'fisher',
    'ranksum'
]

def read_ame(ame_file, method='ranksum'):
    """ Parse the text output of AME into a data frame.

    Please see the documentation for more about AME:
        http://meme-suite.org/doc/ame.html

    N.B. The parsing code here is brittle. If the results seem like nonsense,
    then it is possible the output format of AME has changed, and this code
    needs to be updated.

    The was written based on output from ame version 4.11.2.

    The options to ame which may affect the output were:
        --scoring avg 
        --method ranksum 

    Parameters
    ----------
    ame_file: string
        The path to the txt output of AME

    method: string. see meme_utils.AME_METHODS for allowed methods
        The method used when running AME. 

    Returns
    -------
    ame_df: pd.DataFrame
        A data frame with the following columns:
            motif_name: the name of the motif (from the meme files)
            num_seqs: the number of sequences used in the test
            {left,right,twotailed}_pvalue: the respective p-values
            corrected_{left,right,twotailed}_pvalue: the respective p-values
                after bonferroni correction
            u_value: an "association score" (see the paper for details)
    """
    if method not in AME_METHODS:
        msg = "[meme_utils.read_ame]: invalid method: {}".format(method)
        raise ValueError(msg)

    ame_df = []

    with open(ame_file) as f:
        for line in f:
            if method == 'ranksum':
                res = _parse_ame_ranksum_line(line)
            elif method == 'fisher':
                res = _parse_ame_fisher_line(line)
                        
            ame_df.append(res)
            
    ame_df = collection_utils.remove_nones(ame_df)
    ame_df = pd.DataFrame(ame_df)
    return ame_df