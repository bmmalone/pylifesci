""" This module contains minor helper functions for pyensembl.

Please see the [pyensembl documentation](https://pyensembl.readthedocs.io/en/latest/index.html)
for more details.
"""
import logging
logger = logging.getLogger(__name__)

import numpy as np
import pandas as pd

import pyensembl

import pyllars.logging_utils as logging_utils
import pyllars.collection_utils as collection_utils

###
# Identifier mapping helpers
###
def get_transcript_ids_of_gene_id_df(gene_id:str, ensembl:pyensembl.Genome, raise_on_error:bool=False):
    """ Extract all transcript ids associated with the given gene.

    Parameters
    ----------
    gene_id: string
        The gene identifier

    ensembl: pyensembl.Genome
        The annotations

    raise_on_error: bool
        Whether to raise an error if the gene_id is not present in the
        annotations database

    Returns
    -------
    transcript_gene_id_df: pd.DataFrame
        A dataframe with columns to map between transcripts and genes. Its
        columns are:

            transcript_id
            gene_id
    """
    try:
        transcripts = ensembl.transcript_ids_of_gene_id(gene_id)
    except ValueError as ve:
        msg = ("['pyensembl_utils.get_transcript_ids_of_gene_id_df]: could not "
            "find gene id in database: {}".format(gene_id))
        if raise_on_error:
            raise ValueError(msg) from ve
        else:
            logger.warning(msg)
            return None
    
    ret = [
        {'transcript_id': t, 'gene_id': gene_id}
            for t in transcripts
    ]
    
    return ret

def get_transcript_ids_of_gene_ids(gene_ids, ensembl):
    """ Extract the transcript ids to which the given gene ids match.

    N.B. It is possible for a single transcript id to match to multiple genes.

    Parameters
    ----------
    gene_ids: iterable of strings
        The gene ids

    ensembl: pyensembl.Genome
        The annotations

    Returns
    -------
    transcript_gene_id_df: pd.DataFrame
        A dataframe with columns to map between transcripts and genes. Its
        columns are:

            transcript_id
            gene_id
    """
    id_mapping = parallel.apply_iter_simple(gene_ids, 
        get_transcript_ids_of_gene_id_df, ensembl)
    id_mapping = collection_utils.flatten_lists(id_mapping)
    id_mapping = collection_utils.remove_nones(id_mapping)
    id_mapping = pd.DataFrame(id_mapping)
    return id_mapping

def get_gene_ids_of_transcript_id(transcript_id:str, ensembl:pyensembl.Genome,
        raise_on_error:bool=False):
    """ Extract all gene ids associated with the given transcript.

    Parameters
    ----------
    transcript_id: string
        The transcript identifier

    ensembl: pyensembl.Genome
        The annotations

    raise_on_error: bool
        Whether to raise an exception if the transcript id is not found in the
        annotations database

    Returns
    -------
    transcript_gene_id_df: pd.DataFrame
        A dataframe with columns to map between transcripts and genes. Its
        columns are:

            transcript_id
            gene_id
    """
    try:
        gene_name = ensembl.gene_name_of_transcript_id(transcript_id)
        gene_ids = ensembl.gene_ids_of_gene_name(gene_name)
    except ValueError as ve:
        msg = ("['pyensembl_utils.get_gene_ids_of_transcript_id]: could not "
            "find transcript id in database: {}".format(transcript_id))
        if raise_on_error:
            raise ValueError(msg) from ve
        else:
            logger.warning(msg)
            return None

    
    ret = [
        {'transcript_id': transcript_id, 'gene_id': g} 
            for g in gene_ids
    ]
    
    return ret

def get_gene_ids_of_transcript_ids(transcript_ids, ensembl):
    """ Extract the gene ids to which the given transcript ids match.

    N.B. It is possible for a single transcript id to match to multiple genes.

    Parameters
    ----------
    transcript_ids: iterable of strings
        The transcript ids

    ensembl: pyensembl.Genome
        The annotations

    Returns
    -------
    transcript_gene_id_df: pd.DataFrame
        A dataframe with columns to map between transcripts and genes. Its
        columns are:

            transcript_id
            gene_id
    """    
    id_mapping = parallel.apply_iter_simple(transcript_ids, 
        get_gene_ids_of_transcript_id, ensembl)
    id_mapping = collection_utils.flatten_lists(id_mapping)
    id_mapping = collection_utils.remove_nones(id_mapping)
    id_mapping = pd.DataFrame(id_mapping)
    return id_mapping

def get_gene_name_of_transcript_id(transcript_id:str, ensembl:pyensembl.Genome,
        raise_on_error:bool=False):
    """ Extract the gene name (symbol) for this transcript id.

    The difference between this function and gene_name_of_transcript_id is that
    this function will (optionally) issue a warning rather than raise an
    exception for transcript ids not in the database.

    Parameters
    ----------
    transcript_id: string
        The transcript identifier (e.g., "ENSMUST00000035194")

    ensembl: pyensembl.Genome
        The annotation database

    raise_on_error: bool
        Whether to issue a warning (False) or raise a ValueError (True) if the
        transcript identifier is not in the annotation database

    Returns
    -------
    gene_name: string
        The gene name (also called gene symbol, e.g., "Mapkapk3")

    --- OR ---

    None, if the transcript id is not in the database of annotations
    """
    gene_name = None
    try:
        gene_name = ensembl.gene_name_of_transcript_id(transcript_id)
    except ValueError as ve:
        msg = ("[pyensembl_utils.get_gene_name_of_transcript_id]: could not "
            "find match for transcript id: {}".format(transcript_id))

        if raise_on_error:
            raise ValueError(msg) from ve
        else:
            logger.warning(msg)

    return gene_name

###
# Genome, etc., helpers
###
def get_genome(reference_name, gtf, transcript_fasta=None, logging_args=None,
        annotation_name='ensembl', **kwargs):
    """ Retrieve the pyensembl annotations associated with the given reference.

    The script also creates the database (with the "index" method) if it does
    not already exist.

    This function is largely a wrapper around the pyensembl.Genome constructor,
    so please see it for more details.

    Parameters
    ----------
    reference_name: string
        The identifier for the reference

    gtf: string
        The path to the GTF annotation file

    transcript_fasta: string (optional)
        The path to the fasta file with transcript sequences. The fasta keys
        must match the transcript_id in the GTF file

    logging_args: argparse.Namespace
        pyensembl appears to change several logging levels while opening the
        annotation database (sometimes?). If the logging arguments are given,
        then they will be restored after opening the database.
        
    annotation_name, kwargs:
        Other options to pass to the pyensembl constructor
    """
    ensembl = pyensembl.Genome(
        reference_name=reference_name,
        gtf_path_or_url=gtf,
        transcript_fasta_paths_or_urls=transcript_fasta,
        annotation_name=annotation_name,
        **kwargs
    )
    
    # this will create the database if needed
    ensembl.index()

    if logging_args is not None:
        logging_utils.update_logging(logging_args)

    return ensembl

def concat_protein_sequences(genome=pyensembl.ensembl_grch38, _cache={}):
    """ Create a long string with the longest protein-coding transcript
    version of each gene in the given genome.
    
    https://github.com/openvax/mhc2-data/blob/master/iedb-notebook/Generate%20decoys.ipynb
    """
    if genome in _cache:
        return _cache[genome]
    genes = [gene for gene in genome.genes() if gene.biotype == "protein_coding"]
    all_protein_sequences = []
    for gene in genes:
        curr_transcripts = [
            t for t in gene.transcripts if t.biotype == "protein_coding" 
        ]
        protein_sequences = [
            t.protein_sequence for t in curr_transcripts 
            if  t.protein_sequence is not None 
            and "*" not in t.protein_sequence
        ]
        if len(protein_sequences) == 0:
            continue 
        all_protein_sequences.append(max(protein_sequences, key=len))
    long_str =  "".join(all_protein_sequences)
    _cache[genome] = long_str
    return long_str

def generate_decoys(
        df_epitopes, 
        max_length=None,
        decoys_per_hit=100, 
        genome=pyensembl.ensembl_grch38,
        length_pseudocount=2,
        length_column='epitope_length',
        sequence_column='epitope_sequence'):
    """ Create random, "natural-looking" peptides ("decoys") with length characteristics
    similar to those in the "true" given examples.
    
    https://github.com/openvax/mhc2-data/blob/master/iedb-notebook/Generate%20decoys.ipynb
    """
    
    # calculate the total number of decoys we want
    n_hits = df_epitopes.shape[0]
    n_decoys = n_hits * decoys_per_hit
    
    # create an empirical distribution of lengths
    lengths = df_epitopes[length_column].values
    lengths_set = set(lengths)
    if max_length is None:
        max_length = max(lengths_set)
        
    length_counts = np.ones(max_length, dtype=int) * length_pseudocount
    for n in lengths:
        length_counts[n - 1] += 1
    length_probs = length_counts / length_counts.sum()
    
    # generate a "random" string of all natural sequences
    long_protein_str = concat_protein_sequences(genome)
    
    # randomly sample the starting positions
    max_pos = len(long_protein_str)
    random_indices = np.random.randint(0, max_pos, n_decoys, dtype="int64")
    
    # randomly sample the lengths from the empirical distribution
    random_lengths = np.random.choice(
        a=1 + np.arange(max_length),
        size=n_decoys,
        replace=True,
        p=length_probs)
    
    # collect the random samples
    decoys = []
    for i, n  in zip(random_indices, random_lengths):
        decoys.append(long_protein_str[i:i+n])
        
    df_decoys = pd.DataFrame()
    df_decoys[sequence_column] = decoys
    df_decoys[length_column] = random_lengths
        
    return df_decoys