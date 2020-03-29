#! /usr/bin/env python3

import argparse
import logging
import pandas as pd

import lifesci.mygene_utils as mygene_utils
import pyllars.logging_utils as logging_utils
import pyllars.pandas_utils as pd_utils

logger = logging.getLogger(__name__)

default_filetype = 'AUTO'
filetype_choices = ['AUTO', 'csv', 'excel', 'hdf5']

default_sheet = 0
default_column = None
default_sep = ","

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script creates a report for a list of Entrez or Entrez "
        "gene identifiers. In particular, it extracts information from Swiss-Prot, "
        "TrEMBL, Interpro, PDB, Pfam, PROSITE, the Gene Ontology, and KEGG. It uses "
        "the mygene.info service to collect the annotations.")

    parser.add_argument('filename', help="The name of the file")
    parser.add_argument('out', help="The output file. It will contain the same information "
        "as the input file with the gene annotations appended as new columns. The format "
        "is the same as the input format.")

    parser.add_argument('-f', '--filetype', help="The format of the input file. By default, "
        "the format will be guessed based on the file extension.", choices=filetype_choices,
        default=default_filetype)

    parser.add_argument('--sep', help="The spearator in the file (if csv)", 
        default=default_sep)

    parser.add_argument('-s', '--sheet', help="The name of the sheet (for excel files) "
        "or key (for hdf5 files) from which the gene list is extracted. By default, the "
        "first sheet in an excel file is used. This argument is not used for csv files.",
        default=default_sheet)

    parser.add_argument('-c', '--column', help="The name of the column (or key in hdf5) "
        "from which the gene list is extracted. By default, the first column in the "
        "extracted data frame is used.", default=default_column)

    parser.add_argument('--do-not-compress', help="If this flag is present and the file "
        "is in csv format, then the output will not be compressed. By default, the output "
        "is compressed.", action='store_true')
    
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "Reading the file"
    logger.info(msg)

    df = pd_utils.read_df(args.filename, filetype=args.filetype, sheet=args.sheet, sep=args.sep)

    msg = "Extracting gene identifiers"
    logger.info(msg)

    if args.column is None:
        args.column = df.columns[0]

    gene_ids = df[args.column]

    msg = "Pulling information from mygene.info"
    logger.info(msg)

    res_df = mygene_utils.query_mygene(gene_ids)

    msg = "Joining results to original input"
    logger.info(msg)

    res_df = df.merge(res_df, left_on=args.column, right_on='gene_id', how='inner')

    msg = "Writing output"
    logger.info(msg)

    pd_utils.write_df(res_df, args.out, filetype=args.filetype, sheet=args.sheet,
        do_not_compress=args.do_not_compress, index=False)

    msg = "Finished"
    logger.info(msg)

if __name__ == '__main__':
    main()
