# biopython helpers

import logging
logger = logging.getLogger(__name__)

import numpy as np
import os
import pandas as pd
import re
import shutil
import sys

import pyllars.string_utils as string_utils

logger = logging.getLogger(__name__)

def read_bitseq_alphas(filename, header=None, comment='#', sep=' '):
    """ This function reads the BitSeq m_alphas file into a data frame. In
        particular, the file is assumed to have three columns: theta, alpha
        and beta, as estimated by BitSeq; no header, comment lines beginning
        with '#' and to be space-separated.
            theta, 
        
        Args:
            filename (string): the path to the m_alphas file

            header, comment, sep (strings): these arguments are passed through
                to pandas.read_csv.

        Returns:
            pandas.DataFrame: a data frame with the estimated values as columns
    """
    bitseq_alphas = pd.read_csv(filename.strip(), header=header, comment=comment, sep=sep)
    new_cols = []
    new_cols.append('theta')
    new_cols.append('alpha')
    new_cols.append('beta')
    bitseq_alphas.columns = new_cols

    bitseq_alphas.drop(bitseq_alphas.head(1).index, inplace=True)
    bitseq_alphas.reset_index(drop=True, inplace=True)

    return bitseq_alphas

def read_bitseq_means(filename, names=['mean', 'variance'], comment='#', header=None, sep=' '):
    """ This function reads "mean" files output by getVariance from the BitSeq package.
        This mostly just wraps pd.read_csv; the defaults for the parameters handle the
        standard output format for getVariance, regardless of the quantity estimated
        (RPKM, etc.).

        This data frame should presumably be joined with the dataframe from 
        read_bitseq_tr_file row-wise.

        Args:
            filename (string): the path to the transcript_info file

            names, header, comment, sep (strings): these arguments are passed through
                to pandas.read_csv. N.B. "names" should be a list-like of strings

        Returns:
            pandas.DataFrame: a data frame with the transcript information as 
                columns
    """
    mean_df = pd.read_csv(filename.strip(), comment=comment, names=names, header=header, sep=sep)
    return mean_df

def read_bitseq_tr_file(filename, names=['source_name', 'transcript_id', 
        'length', 'effective_length'], header=None, comment='#', sep=' '):

    """ This function reads the BitSeq transcript_info file into a data frame.
        The file is assumed to contain four columns: source_name, 
        transcript_name, length and effective_length.

        Args:
            filename (string): the path to the transcript_info file

            names, header, comment, sep (strings): these arguments are passed through
                to pandas.read_csv. N.B. "names" should be a list-like of strings

        Returns:
            pandas.DataFrame: a data frame with the transcript information as 
                columns
    """
    transcript_info = pd.read_csv(filename.strip(), names=names, header=header, comment=comment, sep=sep)
    return transcript_info

def read_maxquant_peptides_file(filename, names=None, header='infer', 
        comment='#', sep='\t'):

    """ This function reads the peptides.txt file produced by MaxQuant into a
        data frame. By default, the file is assumed to be tab-delimited, and
        the first row is used as the column names.
    
        Args:
            filename (string): the path to the peptides.txt file

            names, header, comment, sep (strings): these arguments are passed through
                to pandas.read_csv. N.B. "names" should be a list-like of strings or
                None.

        Returns:
            pandas.DataFrame: a data frame with the peptide information
    """
    peptides = pd.read_csv(filename.strip(), names=names, header=header, comment=comment, sep=sep)
    return peptides

def read_protein_digestion_simulator_file(filename, names=None, header='infer',
        comment='#', sep='\t'):

    """ This function reads the output of the Protein Digestion Simulator program
        (https://omics.pnl.gov/software/protein-digestion-simulator). By default,
        the file is assumed to  be tab-delimited and the first row is used as
        column names.

        Args:
            filename (string): the path to the file

            names, header, comment, sep (strings): these arguments are passed through
                to pandas.read_csv. N.B. "names" should be a list-like of strings or
                None.

        Returns:
            pandas.DataFrame: a data frame with the peptide information
    """
    peptides = pd.read_csv(filename.strip(), names=names, header=header, comment=comment, sep=sep)
    return peptides

def get_uniprot_nt_lengths(uniprot_file, min_nt_length=0, types=None):
    """ This function parses a file exported from UniProt giving protein
        information. It then filters based on the 'Status' types, and
        finally filters based on the nucleotide lengths of the sequences.
        It returns a list of lengths.
        
        Args:
            uniprot_file (string) : the path to the file
            
            min_nt_length (int) : the minimum length to include in the result
            
            types (None, or list of strings) : the types to keep after
                filtering; or, if None is given, no filtering will be applied
                
            The available types seem to be:
                * 'unreviewed'
                * 'reviewed'
                * 'unreviewed and UniParc'
                * 'reviewed and UniParc'
                * 'partially reviewed and UniParc'
                * 'partially reviewed'
                
        Returns:
            np.array of ints : the lengths (in nucleotides) of the
                proteins remaining after filtering
    """
    uniprot = pd.read_csv(uniprot_file, sep='\t')

    if types is not None:
        m_types = uniprot['Status'].isin(types)
        uniprot = uniprot[m_types]
        
    uniprot_nt_lengths = uniprot['Length'] * 3

    m_uniprot_nt_lengths = uniprot_nt_lengths > min_nt_length
    uniprot_nt_lengths = uniprot_nt_lengths[m_uniprot_nt_lengths]
    return np.array(uniprot_nt_lengths)

IDMAPPING_SELECTED_NAMES = [
   "UniProtKB-AC",
    "UniProtKB-ID",
    "GeneID (EntrezGene)",
    "RefSeq",
    "GI",
    "PDB",
    "GO",
    "UniRef100",
    "UniRef90",
    "UniRef50",
    "UniParc",
    "PIR",
    "NCBI-taxon",
    "MIM",
    "UniGene",
    "PubMed",
    "EMBL",
    "EMBL-CDS",
    "Ensembl",
    "Ensembl_TRS",
    "Ensembl_PRO",
    "Additional PubMed"
]

def read_idmapping_selected_file(idmapping_selected_file):
    """ This function parses the idmapping_selected.tab file distributed for
        mapping identifiers among different databases.

        The file is available here: 
            ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/idmapping

        and the "README" file (section 2) includes the necessary documentation.
    """
    idmapping_selected = pd.read_csv(idmapping_selected_file, sep='\t', 
        header=None, names=IDMAPPING_SELECTED_NAMES)

    return idmapping_selected

GAF_NAMES = [
    "DB",
    "DB_Object_ID",
    "DB_Object_Symbol",
    "Qualifier",
    "GO_ID",
    "DB:Reference",
    "Evidence Code",
    "With (or) From",
    "Aspect",
    "DB_Object_Name",
    "DB_Object_Synonym",
    "DB_Object_Type",
    "Taxon and Interacting taxon",
    "Date",
    "Assigned_By",
    "Annotation_Extension",
    "Gene_Product_Form_ID"
]

def read_gaf_file(gaf_file):
    """ This function parses the gaf files distributed by the gene ontology 
        into a pandas dataframe.

        A readme file (Section 4, in particular) is here:
            ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/README
    """
    gaf = pd.read_csv(gaf_file, sep='\t', header=None, comment='!', 
        names=GAF_NAMES)

    return gaf

EC2GO_NAMES = [
    'EC', 
    'GO_description', 
    'GO'
]

def read_ec2go_file(ec2go_file):
    """ This function parses the ec2go file distributed by EMBL. The file and
        format information are available here: https://www.ebi.ac.uk/GOA/EC2GO

        N.B. This function assumes the "->" have been converted to symbols
        which do not include ">", since that is used to separate the EC term
        from the GO terms.

        For example, this sed command will fix the symbols.

            sed -i -e 's:->:to:' ec2go.txt
    """
    converters = {
        'EC': str.strip,
        'GO_description': str.strip,
        'GO': str.strip
    }

    ec2go_df = pd.read_csv(ec2go_file, sep='>|;', comment='!', engine='python', 
        header=None, names=EC2GO_NAMES, converters=converters)

    return ec2go_df

NODES_NAMES = [
    "tax_id",
    "parent_tax_id",
    "rank",
    "embl_code",
    "division_id",
    "inherited_div_flag",
    "genetic_code_id",
    "inherited_gc_flag",
    "mitochondrial_genetic_code_id",
    "inherited_mgc_flag",
    "genbank_hidden_flag",
    "hidden_subtree_root_flag",
    "comments"
]

def read_ncbi_taxonomy_nodes_file(nodes_file):
    """ This function parses the nodes.dmp file giving the NCBI taxonomy. It returns it
        as a data frame.
    """
    # pylint: disable=anomalous-backslash-in-string
    nodes_df = pd.read_csv(nodes_file, sep="\t\|\t", header=None, names=NODES_NAMES, engine='python')
    return nodes_df

def get_bowtie2_index_files(base_index_name):
    """ This function returns a list of all of the files necessary for a Bowtie2 index
        that was created with the given base_index_name.

        Args:
            base_index_name (string): the path and base name used to create the bowtie2 index

        Returns:
            list of strings: the paths to all of the files expected for a Bowtie2 index
                based on the provided index_name
    """
    bowtie_extensions = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']
    
    bowtie_files = ['{}{}'.format(base_index_name, bowtie_extension) 
        for bowtie_extension in bowtie_extensions]

    return bowtie_files

###
#   These are helper functions based on identifiers I often use
###
def get_mackowiak_id(object_id, seqname, start, end, strand):
    """ This function creates a "Mackowiak"-style identifier for, e.g., open
        reading frames. The form of the identifier is:

            <object_id>_<seqname>:<start>-<end>:<strand>

        example: ENSMUST00000033123_7:46175479-46179843:-

        start should always be less than end, regardless of the strand; that 
        is, they are bed-style coordinates.

        The ids are "Mackowiak"-style in that they follow the scheme used in
        this paper:
            Mackowiak, S., et al. Extensive identification and analysis of 
                conserved small ORFs in animals. Genome Biology, 2015, 16, 179.

        Args:
            object_id (string): an identifier for the object from which this
                feature was extracted, such as a transcript identifier

            seqname, start, end, strand: the relevant genomic information

        Returns:
            string: the mackowiak-style identifier

    """
    mackowiak_id = "{}_{}:{}-{}:{}".format(object_id, seqname, start, end, strand)
    return mackowiak_id

SPLIT_TOKENS = ["_", ":"]

def parse_mackowiak_id(mackowiak_id):
    """ This function extracts the information out of a "Mackowiak"-style
        identifier. Please see the description for "get_mackowiak_id" to see
        the expected format and information.

        Args:
            mackowiak_id (string): the identifier

        Returns:
            (object_id, seqname, start, end, strand): parsed from the id
    """
    (object_id, seqname, start_end, strand) = string_utils.split(mackowiak_id, SPLIT_TOKENS)

    start, end = start_end.split("-")

    return object_id, seqname, int(start), int(end), strand