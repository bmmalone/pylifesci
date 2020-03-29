#! /usr/bin/env python3

import argparse
import pandas as pd
import re
import yaml

import pyensembl

import misc.utils as utils

import logging
import misc.logging_utils as logging_utils
logger = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Given a meme motif file, extract the gene names and map "
        "them to ensembl identifiers using pyensembl. The pyensembl database "
        "information can be given either in a yaml config file or as command "
        "line options. The yaml config file values have precedence over the "
        "command line options.")

    parser.add_argument('meme', help="The meme file")
    parser.add_argument('out', help="The output file")

    parser.add_argument('-c', '--config', help="The yaml config file. If "
        "given, this should include keys 'genome_name' and 'gtf'. Otherwise, "
        "they may be specified using the respective command line options.",
        default=None)

    parser.add_argument('-n', '--genome-name', help="The genome_parameter for "
        "retrieving the pyensembl database", default=None)
    parser.add_argument('-g', '--gtf', help="The gtf file for pyensembl",
        default=None)
    
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    # if the config file was given, use any values in it to replace those
    # passed on the command line
    if args.config is not None:
        msg = "Reading config file"
        logger.info(msg)
        config = yaml.load(open(args.config))

        args.genome_name = config.get('genome_name', args.genome_name)
        args.gtf = config.get('gtf', args.gtf)

    msg = "genome_name: {}".format(args.genome_name)
    logger.debug(msg)

    msg = "gtf: {}".format(args.gtf)
    logger.debug(msg)

    msg = "Loading pyensembl database"
    logger.info(msg)

    ensembl = pyensembl.Genome(
        reference_name=args.genome_name,
        annotation_name="ensembl", 
        gtf_path_or_url=args.gtf
    )

    # this will create the database if needed
    ensembl.index()

    msg = "Parsing motif gene names"
    logger.info(msg)

    # a line from CISBP looks like:
    #   MOTIF M002_0.6 (Ankhd1)_(Homo_sapiens)_(RBD_1.00)

    all_motifs = []

    motif_re = ("\((?P<gene_name>[^\)]+)\)_\((?P<species>[^\)]+)\)_"
        "\((?P<rbd_score>[^\)]+)\)")
    motif_re = re.compile(motif_re)


    with open(args.meme) as meme_f:
        for line in meme_f:
            if line.startswith("MOTIF"):
                (key, motif_name, info) = line.split()
                
                m = motif_re.match(info)
                
                if m is None:
                    msg = ("Could not parse gene name. Guessing the entire "
                        "string is the gene name: '{}'.".format(info))
                    logger.warning(msg)
                    gene_name = info            
                else:            
                    gene_name = m.group("gene_name")
                
                try:
                    ensembl_ids = ensembl.gene_ids_of_gene_name(gene_name)
                except ValueError:
                    msg = ("Could not find Ensembl identifier for gene_name: "
                        "'{}'".format(gene_name))
                    logger.warning(msg)
                    ensembl_ids = [gene_name]
                
                for ensembl_id in ensembl_ids:
                    motif = {
                        "motif_name": motif_name,
                        "gene_name": gene_name,
                        "ensembl_id": ensembl_id
                    }

                    all_motifs.append(motif)

    msg = "Joining motif gene names into large data frame"
    logger.info(msg)
    all_motifs_df = pd.DataFrame(all_motifs)

    msg = "Writing motifs to disk"
    logger.info(msg)
    utils.write_df(all_motifs_df, args.out, index=False)


if __name__ == '__main__':
    main()
