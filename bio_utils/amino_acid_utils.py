""" This module contains various helpers for working with amino acid sequences.

The `pyllars.string_utils` module also includes many helpers for encoding
amino acid sequences.
"""

import logging
logger = logging.getLogger(__name__)

import io
import joblib
import numpy as np
import pandas as pd
import pathlib
import tqdm

###
# sequence logos
###
from PIL import Image
import weblogo

# openvax
from pepdata.amino_acid_alphabet import canonical_amino_acid_letters as aa_letters

###
# types
###
from typing import Iterable, List, Mapping, Optional, Sequence, Union

path_type = Union[str, pathlib.Path]
aa_encoding_map_type = Mapping[str, np.ndarray]
np_array_or_list = Union[np.ndarray,List]

###
# constans
###
_default_invalid_aa = ['B', 'J', 'O', 'U', 'X', 'Z', '\*']

###
# Loading and creating amino acid property maps
###
def _get_concatenated_properties(
        aa:str,
        aa_properties:Mapping[str, aa_encoding_map_type],
        aa_property_groups:Sequence[str]) -> np.ndarray:
    """ Concatenate the numeric properties for `aa` from all
    `aa_property_groups`
    
    Parameters
    ----------
    aa : str
        The amino acid
        
    aa_properties : typing.Mapping[str, typing.Mapping[str, numpy.ndarray]]
        A map from the property group names to a second mapping for that
        property group. The inner mapping maps from the amino acid to its
        properties for that group.
        
    aa_property_groups : typing.Sequence[str]
        The name of the property groups
        
    Returns
    -------
    concatenated_properties : numpy.ndarray
        The concatenated properties for `aa` from all `aa_property_groups`
    """
    
    concatenated_properties = np.concatenate([
        aa_properties[g][aa]
            for g in aa_property_groups
    ])

    return concatenated_properties

def create_aa_encoding_map(
        aa_properties:Mapping[str, aa_encoding_map_type],
        aa_property_groups:Sequence[str],
        include_j:bool=True,
        include_x:bool=True) -> aa_encoding_map_type:
    """ Create a map from each amino acid to its properties from all
    `aa_property_groups`
    
    Parameters
    ----------
    aa_properties : typing.Mapping[str, typing.Mapping[str, numpy.ndarray]]
        A map from the property group names to a second mapping for that
        property group. The inner mapping maps from an amino acid to its
        numeric properties for that group.
        
    aa_property_groups : typing.Sequence[str]
        The name of the property groups

    include_{j,x} : bool
        Whether to include 'J' ('X') as a "default amino acid" which indicates
        missing values
        
    Returns
    -------
    aa_encoding_map : typing.Mapping[str, numpy.ndarray]
        A map from each amino acid to its concatenated properties from all
        `aa_property_groups`
    """

    aa_keys = aa_letters.copy()

    if include_j:
        aa_keys = aa_keys + ['J']

    if include_x:
        aa_keys = aa_keys + ['X']

    aa_property_groups = sorted(aa_property_groups)    
    aa_property_map = {
        aa: _get_concatenated_properties(aa, aa_properties, aa_property_groups)
            for aa in aa_keys
    }
    
    return aa_property_map


def get_invalid_aa_mask(s_sequences, invalid_aa=_default_invalid_aa):
    """ Find all of the sequences which contain invalid amino acids

    The "standard" non-canonical amino acids are the default for
    `invalid_aa`, that is: B, J, O, U, X, Z, *

    Parameters
    ----------
    s_sequences : pd.Series
        A series (for example, a column from a data frame) containing the amino
        acid sequences.

    invalid_aa : list of chars
        The invalid amino acid characters

    Returns
    -------
    m_invalid : boolean mask
        The sequences which contain invalid amino acid characters are
        marked by a `True`.
    """
    invalid_aa_re = "|".join(invalid_aa)
    m_invalid = s_sequences.str.contains(invalid_aa_re)
    return m_invalid

    
def get_sequences_logo(
        sequences=None,
        probability_matrix=None,
        prior=None,
        alphabet=None,
        formatter='pdf',
        as_pillow=False,
        **kwargs):
    """ Create a logo based on a set of sequences.

    Only one of `sequences` and `probability_matrix` should be given.

    By default, this function assumes that `sequences` would be a list of
    protein sequences, while `probability_matrix` would be for DNA motifs. If
    this is not the case, then the appropriate `alphabet` should be given.

    The most likely values for `alphabet` are:
    * weblogo.seq.unambiguous_dna_alphabet
    * weblogo.seq.unambiguous_rna_alphabet
    * weblogo.seq.unambiguous_protein_alphabet

    Parameters
    ----------
    sequences : iterable of strings
        The sequences

    probability_matrix: 2d np.array
        The proability_matrix of a motif namedtuple

    prior : 1d np.array where len(prior) == len(alphabet)
        A prior to adjust the probabilities as each position

    alphabet : corebio.seq.Alphabet, or None
        The alphabet for creating the motif. If None is given, then the basic
        protein alphabet is used.

    formatter : string
        The format of the image. It must be present in the
        weblogolib.formatters dictionary.

    as_pillow : bool
        Whether to return the raw image data (False, default) or as a
        pillow image (True)

    kwargs : key=value pairs
        Other keyswords to control the image. See weblogolib.LogoOptions for
        complete details. A few likely keywords are:
        - unit_name: for example, "bits" (default) or "probability"
        - fineprint: some text to show at the bottom (default: "WebLogo 3.5.0")
        - show_fineprint: whether to show the fineprint at all
        - yaxis_scale: this appears to control the yaxis maximum
        - stack_width: how wide to make each letter
        - stack_aspect_ratio: the ratio of letter width to height

    Returns
    -------
    logo : binary string or pillow Image
        The raw image data. It can be written to a file opened in binary mode
        ("wb") or opened as a pillow image ("Image.open(io.BytesIO(logo))").

        If `as_pillow` is `True`, then the pillow image is already created
        and returned.
    """
    # ensure we have exactly one of `sequences` and `probability_matrix`
    seq_is_none = (sequences is None)
    pm_is_none = (probability_matrix is None)

    if (seq_is_none and pm_is_none):
        msg = ("[get_sequences_logo] exactly one of `sequences` and "
            "`probability_matrix` must be given")
        raise TypeError(msg) 

    if alphabet is None:
        if sequences is not None:
            alphabet = weblogo.seq.unambiguous_protein_alphabet
        elif probability_matrix is not None:
            alphabet = weblogo.seq.unambiguous_dna_alphabet

    data = None
    if sequences is not None:

        # create the SeqList

        # it seems like there should be a better way to to this...
        sequences = '\n'.join(sequences)
        f = io.StringIO(sequences)
        sequences = weblogo.seq_io.array_io.read(f)
        sequences.alphabet = alphabet
        data = weblogo.LogoData.from_seqs(sequences, prior)
    else: # so we have a probability matrix
        data = weblogo.LogoData.from_counts(
            alphabet, probability_matrix, prior
        )

    options = weblogo.LogoOptions(**kwargs)
    logo_format = weblogo.LogoFormat(data, options)
    formatter = weblogo.formatters[formatter]

    logo = formatter(data, logo_format)

    if as_pillow:
        logo = Image.open(io.BytesIO(logo))

    return logo