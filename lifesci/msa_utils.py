""" Utilities for creating and working with multiple sequence alignments.

Many of these functions are just light wrappers around BioPython alignment utilities:

    https://biopython.org/docs/1.76/api/Bio.Align.Applications.html

The visualization helpers are largely wrappers around the ETE3 toolkit:
    http://etetoolkit.org/docs/latest/tutorial/tutorial_phylogeny.html
"""
import logging

logger = logging.getLogger(__name__)

import pathlib

import Bio.Align.Applications
import ete3


def create_mafft_alignment_tree(
    mafft_unaligned_fasta_filename: pathlib.Path,
    mafft_aligned_fasta_filename: pathlib.Path,
) -> pathlib.Path:
    """Align the given sequences using MAFFT

    **N.B.** This function also creates a phylogenetic tree with the name:
    f"{mafft_unaligned_fasta_filename}.tree". This name is automatically set by MAFFT,
    and it cannot be changed.

    Parameters
    ----------
    mafft_unaligned_fasta_filename : pathlib.Path
        The fasta file containing the input sequences to align

    mafft_aligned_fasta_filename : pathlib.Path
        The output fasta file containing the aligned sequences

    Returns
    -------
    mafft_phylo_tree_filename : pathlib.Path
        The output file with the phylogenetic tree
    """
    mafft_cmd = Bio.Align.Applications.MafftCommandline(
        input=mafft_unaligned_fasta_filename, auto=True, treeout=True
    )

    logger.info("Calling mafft. command: '%s'", mafft_cmd)

    stdout, stderr = mafft_cmd()

    with open(mafft_aligned_fasta_filename, "w") as out:
        out.write(stdout)

    # this is set by mafft
    mafft_phylo_tree_filename = f"{mafft_unaligned_fasta_filename}.tree"
    return mafft_phylo_tree_filename


def draw_msa_phylo_tree(
    mafft_aligned_fasta_filename: pathlib.Path,
    mafft_phylo_tree_filename: pathlib.Path,
    phylo_tree_out: pathlib.Path,
) -> ete3.PhyloTree:
    """Draw the MSA and phylogenetic tree based on the given MAFFT alignments

    **N.B.** The species/identifiers in the alignment fasta and tree files must match
    exactly. MAFFT replaces may "odd" characters (such as `.`) with underscores. If
    the identifiers do not match across the two files, then a warning similar to the
    following will be printed, and the image will not include the MSA.

    > terminal nodes could not be found in the alignment.

    Parameters
    ----------
    mafft_{aligned_fasta,phylo_tree}_filename : pathlib.Path
        The respective output files from MAFFT

    phylo_tree_out : pathlib.Path
        The path to the output image

    Returns
    -------
    tree : ete3.PhyloTree
        The phylogenetic tree object
    """
    newick = ""  # newick is the name of the phylo tree format
    with open(mafft_phylo_tree_filename, "r") as f:
        lines = f.readlines()
        for l in lines:
            l = l.strip()

            # MAFFT prepends a number and underscore to each species/identifier in its
            # tree output file. This is not compatible with ECE3, so remove those.
            if l[0].isdigit():
                l = l.split("_", maxsplit=1)
                l = l[1]

            newick += l

    tree = ete3.PhyloTree(newick)
    tree.link_to_alignment(alignment=mafft_aligned_fasta_filename, alg_format="fasta")
    tree.render(phylo_tree_out)

    return tree
