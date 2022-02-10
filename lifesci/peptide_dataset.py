""" This module contains a class for loading peptide datasets, including
standard filtering operations.
"""
import logging
_logger = logging.getLogger(__name__)

import numpy as np
import pandas as pd
import pyllars.pandas_utils as pd_utils

import lifesci.amino_acid_utils as amino_acid_utils

from typing import Mapping, Optional, Sequence, Set

_DEFAULT_NAME = "PeptideDataset"

def _length_filter(valid_lengths:Set[int], peptides:pd.Series) -> np.ndarray:
    m_len = peptides.str.len().isin(valid_lengths)
    return m_len

def _9mers_only_filter(peptides):
    m_len = _length_filter({9}, peptides)
    return m_len
    
def _15mers_only_filter(peptides):
    m_len = _length_filter({15}, peptides)
    return m_len

def _standard_aa_only_filter(peptides):
    m_invalid_aa = amino_acid_utils.get_invalid_aa_mask(peptides)
    m_valid_aa = ~m_invalid_aa
    return m_valid_aa

FILTER_MAP = {
    "9mers_only": _9mers_only_filter,
    "15mers_only": _15mers_only_filter,
    "standard_aa_only": _standard_aa_only_filter
}

class PeptideDataset(object):

    def __init__(self,
            filename:str,
            peptide_column:str,
            loading_parameters:Mapping=dict(),
            filters:Sequence[str]=list(),
            name:str=_DEFAULT_NAME):
        self.filename = filename
        self.peptide_column = peptide_column
        self.loading_parameters = loading_parameters
        self.filters = filters
        self.name = name

    def log(self,
            msg:str,
            logger:Optional[logging.Logger]=None,
            level:int=logging.INFO) -> None:
        if logger is None:
            logger = _logger

        msg = "[{}]: {}".format(self.name, msg)
        logger.log(level, msg)

    def load_dataset(self):
        df_peptides = pd_utils.read_df(self.filename, **self.loading_parameters)

        peptides = df_peptides[self.peptide_column]

        if self.filters is not None:
            all_filters = [
                FILTER_MAP[f](peptides) for f in self.filters
            ]

            m_filters = pd_utils.intersect_masks(all_filters)

            df_peptides = df_peptides[m_filters].reset_index(drop=True)

        return df_peptides

    @classmethod
    def load(clazz,
            filename:str,
            peptide_column:str,
            filters:Optional[Sequence[str]]=None,
            loading_parameters:Mapping=dict()) -> pd.DataFrame:
        
        d = clazz(
            filename=filename,
            peptide_column=peptide_column,
            loading_parameters=loading_parameters,
            filters=filters
        )

        df_peptides = d.load_dataset()
        return df_peptides