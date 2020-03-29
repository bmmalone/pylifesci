""" This module contains helpers for working with the mhcnames package.

mhcnames is part of OpenVax, and it is distributed using the Apache 2.0
license. Please see the github page for more details: https://github.com/openvax/mhcnames
"""
import logging
logger = logging.getLogger(__name__)

import numpy as np
import pandas as pd

import mhcnames
import pyllars.collection_utils as collection_utils

from typing import Iterable, List, Optional, Sequence

def get_drb(s:str) -> str:
    """ Parse any "DRB" allele from `s`
    """
    s = s.split("-")    
    drb = None
    
    if s[0].startswith("DRB"):
        drb = s[0]
    elif s[1].startswith("DRB"):
        drb = s[1]
        
    return drb

def get_normalized_allele_name(raw_allele_name:str) -> Optional[str]:
    """ Try to parse `raw_allele_name` to a normalized compact allele name
    
    If `raw_allele_name` cannot be parsed, a warning is logged and `None` is
    returned.
    
    Parameters
    ----------
    raw_allele_name : str
        The unnormalized allele name
        
    Returns
    -------
    compact_allele_name : typing.Optional[str]
        The normalized, compact allele name, or `None`, if the allele name could
        not be parsed.
    """
    ret = None
    
    try:
        ret = mhcnames.compact_allele_name(raw_allele_name)
    except mhcnames.AlleleParseError as ape:
        logger.warning(ape)
        
    return ret

def get_normalized_allele_names(raw_allele_names:Sequence[str]) -> List[str]:
    """ Get a list of normalized allele names from `raw_allele_names`
    
    Parameters
    ----------
    raw_allele_names : typing.Sequence[str]
        The allele names. Please see the `mhcnames` documentation for how different
        forms of allele names are handled.
        
    Returns
    -------
    normalized_allele_names : typing.List[str]
        The normalized allele name for each allele in `raw_allele_names`. Any allele
        names which could not be converted to a normalized form will have `None` as
        its value in this list.
    """
    m_missing = pd.isnull(raw_allele_names)
    unique_allele_names = set(raw_allele_names[~m_missing])
    
    normalized_allele_name_map = {
        a: get_normalized_allele_name(a)
            for a in unique_allele_names
    }
    
    normalized_allele_names = [
        normalized_allele_name_map.get(a, None)
            for a in raw_allele_names
    ]
    
    return normalized_allele_names


def _normalize_allele_name(parsed_allele):
    """ A helper for `_normalize_ab_allele_names` to account for allele families
    """    
    normalized_allele_name = None
    if len(parsed_allele.allele_family) > 0:
        
        normalized_allele_name = "{}{}{}".format(
            parsed_allele.gene,
            parsed_allele.allele_family,
            parsed_allele.allele_code
        )
    else:
        # mice don't have allele families
        normalized_allele_name = "{}{}".format(
            parsed_allele.gene,
            parsed_allele.allele_code
        )
            
    return normalized_allele_name
    

def _normalize_ab_allele_names(ab_allele_names):
    """ A helper for `get_compact_ab_allele_names`
    """
    
    try:
        res = mhcnames.parse_classi_or_classii_allele_name(ab_allele_names)
    except Exception as e:
        msg = str(e)
        logger.warning(msg)
        return None
    
    if len(res) != 2:
        msg = ("Could not parse alpha-beta allele name: {}".format(
            ab_allele_names))
        logger.warning(msg)
        return None
        
    alpha, beta = res
    
    alpha_normalized = _normalize_allele_name(alpha)
    beta_normalized = _normalize_allele_name(beta)
    
    ab_normalized = "{}-{}".format(alpha_normalized, beta_normalized)
    return ab_normalized
    
def get_normalized_ab_allele_names(ab_allele_names:Iterable[str]) -> List[str]:
    """ Get a list of normalized allele names based on the given list
    
    In contract to `get_normalized_allele_names`, this function determines both the
    alpha and beta chain names. Please see the `mhcnames.parse_classi_or_classii_allele_name`
    documentation for chain inference when both chains are not specified.
    
    This function (eventually) uses the `mhcnames.compact_allele_name`
    function to determine the normalized name for each allele.
    
    Parameters
    ----------
    ab_allele_names : typing.Iterable[str]
        The allele names. Please see the `mhcnames` documentation for how different
        forms of allele names are handled.
        
        The only restriction on this input is that `np.unique` should return a list
        of the unique values.
        
    Returns
    -------
    normalized_ab_allele_names : typing.List[str]
        The normalized alpha-beta allele name for each allele in `allele_names`
    """
    
    normalized_ab_allele_name_map = {
        a: _normalize_ab_allele_names(a)
            for a in np.unique(ab_allele_names)
    }
    
    normalized_ab_allele_names = [
        normalized_ab_allele_name_map[a] for a in ab_allele_names
    ]
    
    return normalized_ab_allele_names