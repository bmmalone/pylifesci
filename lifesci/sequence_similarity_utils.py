"""
Utilities for calling and parsing calls to various Blast+ programs, BLOSUM
similarity calclations, and calculating SNebula scores.

See the documentation for BLAST:
    https://blast.ncbi.nlm.nih.gov/Blast.cgi

The blast helpers are largely wrappers around the blast modules from BioPython.
    http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc96
"""
import logging

logger = logging.getLogger(__name__)

import numpy as np
import pandas as pd
import pathlib
import shutil
import tempfile

import pyllars.pandas_utils as pd_utils

import lifesci.amino_acid_utils as aa_utils
import lifesci.fastx_utils as fastx_utils

import Bio
import Bio.pairwise2
from Bio.Blast.Applications import NcbiblastpCommandline

from typing import Optional

_DEFAULT_BLAST_OUT_FORMAT = '"6 qseqid qseq sseqid sseq qcovs pident mismatch evalue bitscore gapopen gaps qstart qend sstart send"'

_DEFAULT_BLAST_OUT_COLUMNS = [
    "query_ID",
    "query_seq_aligned",
    "database_ID",
    "database_seq_aligned",
    "percent_query_coverage",
    "percent_alignement",
    "num_mismatches",
    "e_value",
    "bitscore",
    "gapopen",
    "gaps",
    "query_seq_alignment_start",
    "query_seq_alignment_end",
    "database_seq_alignment_start",
    "database_seq_alignment_end",
]

# N.B. This currently raises a deprecation error, but using `Bio.Align.substitution_matrices`
# is significantly slower than using Bio.SubsMat.MatrixInfo
# import Bio.Align.substitution_matrices
# DEFAULT_SCORING_MATRIX = Bio.Align.substitution_matrices.load("BLOSUM62")
import Bio.SubsMat.MatrixInfo

DEFAULT_SCORING_MATRIX = Bio.SubsMat.MatrixInfo.blosum62

# some very large threshold so we always find some matches
DEFAULT_E_THRESHOLD = 1e10

###
# SNebula-related helpers
###
def _calc_l_xi(l_x, l_i):
    """Calculate the length normalization term for SNebula"""
    num = np.abs(l_i - l_x)
    den = l_x

    l_xi = num / den
    return l_xi


def get_snebula_score(p_x, p_i, scoring_matrix=DEFAULT_SCORING_MATRIX):
    """Calculate the similarity between `p_x` and `p_i` using the SNebula
    formulas

    Parameters
    ----------
    p_{x,i} : str
        The peptide sequences

    scoring_matrix : dict
        A dictionary which maps two-tuples to scores. Each element in the two
        tuples are amino acids, and the value gives the score of aligning the
        first element with the second element.

        Presumably, this is one of the options from `Bio.SubsMat.MatrixInfo`.

    Returns
    -------
    score : float
        The snebula similarity score
    """

    p_xi = Bio.pairwise2.align.globaldx(p_x, p_i, scoring_matrix, score_only=True)

    p_xx = Bio.pairwise2.align.globaldx(p_x, p_x, scoring_matrix, score_only=True)

    l_xi = _calc_l_xi(len(p_x), len(p_i))
    e_term = np.exp(-1 * l_xi)

    val = p_xi / p_xx * e_term
    return val


def _get_snebula_score_helper(row, scoring_matrix=DEFAULT_SCORING_MATRIX):
    mutated_epitope = row["query_seq_aligned_ungapped"]
    reference_epitope = row["database_seq_aligned_ungapped"]

    score = get_snebula_score(mutated_epitope, reference_epitope, scoring_matrix)

    ret = {
        "query_seq_aligned": row["query_seq_aligned"],
        "database_seq_aligned": row["database_seq_aligned"],
        "query_seq_aligned_ungapped": mutated_epitope,
        "database_seq_aligned_ungapped": reference_epitope,
        "database_ID": row["database_ID"],
        "snebula_score": score,
        "num_mismatches": row["num_mismatches"],
    }

    return ret


###
# BLOSUM-related helpers
###
def get_normalized_similarity(
    mutated_epitope, reference_epitope, scoring_matrix=DEFAULT_SCORING_MATRIX
):
    """Find the alignment similarity between the two epitopes

    This function normalizes the score between the two epitopes using the
    similarity score from `mutated_epitope` to itself. This is, the normalized
    similarity is given as: \frac{sim(mut, ref)}{sim(mut, mut)}.

    Please see `get_pairwise_normalized_similarity` for a version which
    normalizes by both epitope sequences.

    Parameters
    -----------
    {mutated,reference}_epitope : str
        The respective epitope sequences. They can have different lengths.

    scoring_matrix : dict
        A dictionary which maps two-tuples to scores. Each element in the two
        tuples are amino acids, and the value gives the score of aligning the
        first element with the second element.

        Presumably, this is one of the options from `Bio.SubsMat.MatrixInfo`.

    Returns
    -------
    score : float
        The normalized score as described above.
    """

    score_reference_mutant = Bio.pairwise2.align.globaldx(
        mutated_epitope, reference_epitope, scoring_matrix, score_only=True
    )

    score_mutant_mutant = Bio.pairwise2.align.globaldx(
        mutated_epitope, mutated_epitope, scoring_matrix, score_only=True
    )

    return score_reference_mutant / score_mutant_mutant


def get_pairwise_normalized_similarity(
    epitope_1, epitope_2, scoring_matrix=DEFAULT_SCORING_MATRIX
):
    """Find the alignment similarity between the two epitopes

    This function normalizes the score between the two epitopes using the
    "self-similarity" of both epitope sequences. This is, the normalized
    similarity is given as: \frac{2*sim(e1, e2)}{sim(e1, e1) + sim(e2, e2)}.

    Please see `get_normalized_similarity` for a version which
    normalizes only by the self-similarity of `epitope_1`.

    Parameters
    -----------
    epitope_{1,2} : str
        The respective epitope sequences. They can have different lengths.

    scoring_matrix : dict
        A dictionary which maps two-tuples to scores. Each element in the two
        tuples are amino acids, and the value gives the score of aligning the
        first element with the second element.

        Presumably, this is one of the options from `Bio.SubsMat.MatrixInfo`.

    Returns
    -------
    score : float
        The normalized score as described above.
    """

    e1_e2_score = Bio.pairwise2.align.globaldx(
        epitope_1, epitope_2, scoring_matrix, score_only=True
    )

    e1_e1_score = Bio.pairwise2.align.globaldx(
        epitope_1, epitope_1, scoring_matrix, score_only=True
    )

    e2_e2_score = Bio.pairwise2.align.globaldx(
        epitope_2, epitope_2, scoring_matrix, score_only=True
    )

    sim = 2 * e1_e2_score / (e1_e1_score + e2_e2_score)
    return sim


###
# BLAST-related helpers
###


def _get_query_db_length_matches(df, length):
    m_query = df["query_seq_len"] == length
    m_db = df["database_seq_len"] == length
    m_match = m_query & m_db
    return m_match


def _get_similarity_helper(row, scoring_matrix=DEFAULT_SCORING_MATRIX):
    mutated_epitope = row["query_seq_aligned_ungapped"]
    reference_epitope = row["database_seq_aligned_ungapped"]

    score = get_normalized_similarity(
        mutated_epitope, reference_epitope, scoring_matrix
    )

    ret = {
        "query_seq_aligned": row["query_seq_aligned"],
        "database_seq_aligned": row["database_seq_aligned"],
        "query_seq_aligned_ungapped": mutated_epitope,
        "database_seq_aligned_ungapped": reference_epitope,
        "database_ID": row["database_ID"],
        "normalized_similarity": score,
        "num_mismatches": row["num_mismatches"],
    }

    return ret


def get_similarity_from_blast_hits(
    sequences,
    blastdb_path: pathlib.Path,
    num_threads: int = 1,
    identifiers=None,
    e_value: float = DEFAULT_E_THRESHOLD,
    scoring_matrix=DEFAULT_SCORING_MATRIX,
    query_output_folder: Optional[pathlib.Path] = None,
    delete_query_output_folder: bool = True,
    discard_duplicate_matches: bool = True,
    ungapped: bool = True,
    qcov_hsp_perc: int = 100,
    include_complete_blast_output: bool = False,
    progress_bar: bool = True,
    gapextend: int = 1,
    gapopen: int = 9,
    word_size: int = 3,
    max_target_seqs: int = 500,
):

    # for logging
    caller = "get_similarity_from_blast_hits"

    if identifiers is None:
        msg = "[{}]: creating identifiers".format(caller)
        logger.info(msg)
        identifiers = ["seq_{}".format(i) for i in range(len(sequences))]

    if query_output_folder is None:
        query_output_folder = tempfile.mkdtemp()

    msg = "[{}]: using temp folder: {}".format(caller, query_output_folder)
    logger.info(msg)
    query_output_folder = pathlib.Path(query_output_folder)

    msg = "[{}]: writing the sequences and identifiers to disk".format(caller)
    logger.info(msg)
    query_filename = query_output_folder / "query.fa"

    seq_ids = list(zip(identifiers, sequences))
    fastx_utils.write_fasta(seq_ids, query_filename, compress=False)

    # create our output path
    blastp_filename = query_output_folder / "blastp.tab"

    blp = NcbiblastpCommandline(
        query=str(query_filename),
        out=str(blastp_filename),
        db=blastdb_path,
        outfmt=_DEFAULT_BLAST_OUT_FORMAT,
        evalue=e_value,
        ungapped=ungapped,
        comp_based_stats="F",
        num_threads=num_threads,
        qcov_hsp_perc=qcov_hsp_perc,
        gapextend=gapextend,
        gapopen=gapopen,
        word_size=word_size,
        max_target_seqs=max_target_seqs,
    )

    msg = "calling blastp. command: '{}'".format(blp)
    logger.info(msg)

    stdout, stderr = blp()

    msg = "[{}]: loading blastp output".format(caller)
    logger.info(msg)
    df_blastp = pd.read_csv(
        blastp_filename, sep="\t", header=None, names=_DEFAULT_BLAST_OUT_COLUMNS
    )

    m_duplicates = np.array([False] * len(df_blastp))
    if discard_duplicate_matches:
        msg = "[{}]: filtering duplicate matches".format(caller)
        logger.info(msg)

        cols = ["query_seq", "database_seq"]
        m_duplicates = df_blastp.duplicated(subset=cols)

    msg = "[{}]: filtering sequences with non-standard aa".format(caller)
    logger.info(msg)
    m_invalid_aa = aa_utils.get_invalid_aa_mask(df_blastp["database_seq_aligned"])

    # keep only matches which pass all filters
    m_to_keep = ~m_duplicates & ~m_invalid_aa
    df_blastp = df_blastp[m_to_keep]
    df_blastp = df_blastp.reset_index(drop=True)

    msg = "[{}]: calculating sequence similarities".format(caller)
    logger.info(msg)

    df_blastp["query_seq_aligned_ungapped"] = df_blastp[
        "query_seq_aligned"
    ].str.replace("-", "")
    df_blastp["database_seq_aligned_ungapped"] = df_blastp[
        "database_seq_aligned"
    ].str.replace("-", "")

    similarities = pd_utils.apply(
        df_blastp,
        # _get_similarity_helper,
        _get_snebula_score_helper,
        scoring_matrix=scoring_matrix,
        progress_bar=progress_bar,
    )

    msg = "[{}]: building result data frame".format(caller)
    logger.info(msg)
    df_similarities = pd.DataFrame(similarities)

    if include_complete_blast_output:
        # the "concat_cols" depends on the exact score used
        # concat_cols = ["normalized_similarity"]
        concat_cols = ["snebula_score"]
        df_similarities = pd.concat([df_blastp, df_similarities[concat_cols]], axis=1)
        df_similarities["query_seq_aligned_len"] = df_similarities[
            "query_seq_aligned"
        ].str.len()
        df_similarities["database_seq_aligned_len"] = df_similarities[
            "database_seq_aligned"
        ].str.len()

    if delete_query_output_folder:
        shutil.rmtree(query_output_folder)

    return df_similarities
