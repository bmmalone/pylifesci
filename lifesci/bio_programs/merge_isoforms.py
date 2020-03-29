#! /usr/bin/env python3

import argparse
import csv
import logging
import os
import re
import numpy as np
import pandas as pd
import shlex
import sys

import lifesci.bio as bio
import lifesci.bed_utils as bed_utils
import lifesci.gtf_utils as gtf_utils
import pyllars.logging_utils as logging_utils
import misc.parallel as parallel
import misc.slurm as slurm

logger = logging.getLogger(__name__)

default_num_groups = 500

default_feature_type = 'CDS'
default_group_attribute = 'gene_id'
default_id_format_str = "{}.merged"

default_chr_name_file = None

def parse_attributes_group(rows):
    res = parallel.apply_df_simple(rows, gtf_utils.parse_gtf_attributes)
    res = pd.DataFrame(res)
    return res

def merge_gene_group(g, transcript_id_field, transcript_id_format_str):

    interval_starts = np.array(g['start'])

    # we must add 1 to the ends because gtf is !closed! at the end, but
    # the merge function expects the end to be open
    interval_ends = np.array(g['end']) + 1

    res = bed_utils.merge_intervals(interval_starts, interval_ends)
    merged_starts, merged_ends, _ = res

    # subtract 1 back off the ends
    merged_ends = merged_ends - 1
    
    # convert back to gtf
    merged_df = pd.DataFrame()
    merged_df['start'] = merged_starts
    merged_df['end'] = merged_ends

    merged_df['seqname'] = g['seqname'].iloc[0]
    merged_df['source'] = g['source'].iloc[0]
    merged_df['feature'] = g['feature'].iloc[0]
    merged_df['score'] = "0"
    merged_df['strand'] = g['strand'].iloc[0]
    merged_df['frame'] = "."
    
    gene_id = None
    if 'gene_id' in g:
        gene_id = g['gene_id'].iloc[0]

    gene_name = None
    if 'gene_name' in g:
        gene_name = g['gene_name'].iloc[0]

    transcript_id = g[transcript_id_field].iloc[0]
    transcript_id = transcript_id_format_str.format(transcript_id)

    attributes = [
        "transcript_id " +
        '"' +
        transcript_id +
        '"'
    ]

    if gene_id is not None:
        attributes.extend([
            "gene_id " +
            '"' +
            gene_id +
            '"'
        ])

    if gene_name is not None:
        attributes.extend([
            "gene_name" +
            '"' +
            gene_name +
            '"'
        ])

    attributes = '; '.join(attributes)    
    merged_df['attributes'] = attributes
    
    return merged_df[gtf_utils.gtf_field_names]

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script merges either the exons or CDS regions of all transcript "
        "isoforms into a single \"super gene isoform\". It does this based on the given "
        "GTF feature type and attribute (with defaults \"CDS\" and \"gene_id\", respectively).")
    parser.add_argument('gtf', help="The GTF file")
    parser.add_argument('out', help="The output (merged) GTF file")

    parser.add_argument('--feature-type', help="The type of features to merge", 
        default=default_feature_type)
    parser.add_argument('--group-attribute', help="The attribute by which the features "
        "will be merged", default=default_group_attribute)

    parser.add_argument('--id-format-str', help="The python format string to "
        "use for creating the \"transcript\" identifiers", 
        default=default_id_format_str)

    parser.add_argument('--chr-name-file', help="If this file is specified, it will "
        "be used to determine the seqname sort order. This should be the "
        "\"chrName.txt\" file created by STAR. If not present, the transcripts "
        "will be sorted alphabetically (1, 10, 11, 2, ..., KL568162.1, MT, X, Y).",
        default=default_chr_name_file)

    parser.add_argument('--add-exons', help="If this flag is given, then all features will "
        "be duplicated, but with the feature type \"exon\". Presumably, this should be given "
        "when \"CDS\" features are merged, and the resulting GTF file will be used by STAR "
        "(or anything else expecting \"exon\"s).", action='store_true')
    
    parser.add_argument('-g', '--num-groups', help="The number of groups into which to split "
        "the features. More groups means the progress bar is updated more frequently but incurs "
        "more overhead because of the parallel calls.", type=int, default=default_num_groups)
    
    slurm.add_sbatch_options(parser)
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    if args.use_slurm:
        cmd = ' '.join(sys.argv)
        slurm.check_sbatch(cmd, args=args)
        return

    msg = "Reading GTF file"
    logger.info(msg)

    gtf_df = gtf_utils.read_gtf(args.gtf)

    msg = "Extracting desired features"
    logger.info(msg)
    m_feature_type = gtf_df['feature'] == args.feature_type
    gtf_feature_df = gtf_df[m_feature_type]

    msg = "Parsing GTF attributes"
    logger.info(msg)

    attributes = parallel.apply_parallel_split(
        gtf_feature_df,
        args.num_cpus, 
        parse_attributes_group,
        progress_bar=True,
        num_groups=args.num_groups
    )

    attributes_df = pd.concat(attributes)
    attributes_df['end'] = attributes_df['end'].astype(int)
    attributes_df['start'] = attributes_df['start'].astype(int)

    msg = "Merging isoforms"
    logger.info(msg)

    gene_features = attributes_df.groupby(args.group_attribute)
    merged_genes = parallel.apply_parallel_groups(
        gene_features,
        args.num_cpus, 
        merge_gene_group,
        args.group_attribute,
        args.id_format_str,
        progress_bar=True
    )

    merged_genes_df = pd.concat(merged_genes)

    if args.add_exons:
        merged_exons = merged_genes_df.copy()
        merged_exons['feature'] = 'exon'
        merged_genes_df = pd.concat([merged_exons, merged_genes_df])

    merged_genes_df['start'] = merged_genes_df['start'].astype(int)

    # now, sort the merged isoforms

    # this is a bit of a hack, because it is actually using the sorting routine
    # for bed data frames

    # we need a dummy 'id' column for sorting, so just use the attributes
    merged_genes_df['id'] = merged_genes_df['attributes']
    merged_genes_df = bed_utils.sort(
        merged_genes_df,
        seqname_order=args.chr_name_file
    )

    # last, drop duplicate rows
    fields = ['seqname', 'source', 'feature', 'start', 'end', 'strand']
    merged_genes_df = merged_genes_df.drop_duplicates(subset=fields)

    gtf_utils.write_gtf(merged_genes_df, args.out, compress=False)

if __name__ == '__main__':
    main()
