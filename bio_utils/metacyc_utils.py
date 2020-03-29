""" Utilities for working with the metacyc database
"""
import logging
logger = logging.getLogger(__name__)

import pandas as pd

import pyllars.validation_utils as validation_utils


import sexpdata
#sexpdata (for parsing the S-exp data: https://github.com/tkf/sexpdata

def parse_metacyc_attribute(line):
    """ This function parses the key and value from a metacyc dat file.
        The format is described here: http://bioinformatics.ai.sri.com/ptools/flatfile-format.html
    """
    key, val = line.strip().split(" - ")
    return key, val

def parse_pathways_species(pathways_file):
    """ This function parses a metacyc pathways.dat file to produce a data frame
        which includes (pathway,species) mappings.

        Here, the species is typically an NCBI taxonomy identifier; in these
        cases, the form is "TAX-<tax id>". Other times, it is a higher-order
        object, such as an order. These id's have the form "ORG-<id>".

        Args:
            pathways_file (string): path to the pathways.dat file

        Returns:
            pd.DataFrame with columns: pathway, species
    """
    pathways_species = []

    cur_pathway = None

    # there are sometimes weird encoding issues here; just skip them
    with open(pathways_file, encoding='utf-8', errors='ignore') as pf:
        for line in pf:
            if line.startswith("UNIQUE-ID"):
                key, val = parse_metacyc_attribute(line)
                cur_pathway = val
                
            elif line.startswith("SPECIES"):
                if cur_pathway is None:
                    continue
                    
                key, val = parse_metacyc_attribute(line)
                pathway_species = {
                    "pathway": cur_pathway,
                    "species": val
                }
                
                pathways_species.append(pathway_species)
                
    pathways_species_df = pd.DataFrame(pathways_species)
    return pathways_species_df

def parse_pathways_ecs(reactions_file):
    """ This function parses a metacyc reactions.dat file to produce a data frame
        which include (pathway,ec) mappings.

        The reactions.dat file assigns an ec-number to each reaction; it also
        gives the pathways to which each reaction belongs. This function joins
        ec numbers and pathways on the reactions to find the mapping.

        Args:
            reactions_file (string): path to the reactions.dat file

        Returns:
            pd.DataFrame with columns: pathway, ec
    """
    pathways_ecs = []

    cur_reaction = None
    cur_ec_number = None
    with open(reactions_file, encoding='utf-8', errors='ignore') as rf:
        for line in rf:
            if line.startswith("UNIQUE-ID"):
                key, val = parse_metacyc_attribute(line)
                cur_reaction = val            
                cur_ec_number = None
                
            elif line.startswith("EC-NUMBER"):
                if cur_reaction is None:
                    continue
                    
                key, val = parse_metacyc_attribute(line)
                cur_ec_number = val
                
            elif line.startswith("IN-PATHWAY"):
                if cur_ec_number is None:
                        continue
                        
                key, val = parse_metacyc_attribute(line)
                pathway_ec = {
                    "ec": cur_ec_number,
                    "pathway": val
                }
                
                pathways_ecs.append(pathway_ec)
                
    pathways_ecs_df = pd.DataFrame(pathways_ecs)
    return pathways_ecs_df

def parse_uniprot_ecs(uniprot_seq_ids_file, raise_on_error=True, logger=logger):
    """ This function parses a metacyc uniprot-seq-ids.dat file to produce a
        data frame which includes (ec,uniprot_id) mappings.

        The uniprot-seq-ids.dat file directly gives the mapping between these
        in lisp S-exp format.

        Args:
            uniprot_seq_ids_file (string): path to the uniprot-seq-ids.dat file

        Returns:
            pd.DataFrame with columns: ec, uniprot_id
    """
    data = sexpdata.load(open(uniprot_seq_ids_file))
    ec_uniprot_list = []

    for entry in data:
        ec = entry[1]
        # we sometimes have lists of size 1 of ECs. It is not clear why this
        # happens, so just parse it out
        if validation_utils.validate_is_sequence(ec, raise_on_invalid=False):
            # if the size is not 1, then this is unexpected
            if len(ec) > 1:
                msg = ("Unexpected sequence when parsing S-exp EC from uniprot_"
                    "seqids_file. {}".format(data))
                if raise_on_error:
                    raise ValueError(msg)
                else:
                    logger.warning(msg)
            ec = ec[0]

        # the EC is sometimes as "Symbol"; this is okay
        if isinstance(ec, sexpdata.Symbol):
            ec = ec.value()
        
        for uniprot_id in entry[2:]:
            # we occasionally parse out lists of size 1 of symbols
            # it is unclear to me when this happens
            if validation_utils.validate_is_sequence(uniprot_id, raise_on_invalid=False):
                # if the size of the list is not 1, then this is unexpected
                if len(uniprot_id) > 1:
                    msg = ("Unexpected sequence when parsing S-exp uniprot_"
                        "seq_ids_file. {}".format(data))
                    if raise_on_error:
                        raise ValueError(msg)
                    else:
                        logger.warning(msg)

                uniprot_id = uniprot_id[0]

            if isinstance(uniprot_id, sexpdata.Symbol):
                uniprot_id = uniprot_id.value()

            ec_uniprot = {
                "ec": ec,
                "uniprot_id": uniprot_id
            }
            
            ec_uniprot_list.append(ec_uniprot)
            
    ec_uniprot_df = pd.DataFrame(ec_uniprot_list)

    return ec_uniprot_df

