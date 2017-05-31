###
#   This module contains minor helper functions for pyensembl.
###
import logging
import pyensembl

logger = logging.getLogger(__name__)

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
    import misc.parallel as parallel
    import misc.utils as utils
    import pandas as pd

    id_mapping = parallel.apply_iter_simple(gene_ids, 
        get_transcript_ids_of_gene_id_df, ensembl)
    id_mapping = utils.flatten_lists(id_mapping)
    id_mapping = utils.remove_nones(id_mapping)
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
    import misc.parallel as parallel
    import misc.utils as utils
    import pandas as pd
    
    id_mapping = parallel.apply_iter_simple(transcript_ids, 
        get_gene_ids_of_transcript_id, ensembl)
    id_mapping = utils.flatten_lists(id_mapping)
    id_mapping = utils.remove_nones(id_mapping)
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
    import misc.logging_utils as logging_utils
    import pyensembl

    ensembl = pyensembl.Genome(
        reference_name=reference_name,
        gtf_path_or_url=gtf,
        transcript_fasta_path_or_url=transcript_fasta,
        annotation_name=annotation_name,
        **kwargs
    )
    
    # this will create the database if needed
    ensembl.index()

    if logging_args is not None:
        logging_utils.update_logging(logging_args)

    return ensembl
