#! /usr/bin/env python3

import matplotlib
matplotlib.use('agg')

import argparse
import matplotlib.ticker as mtick
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns; sns.set(style='white')

import misc.mpl_utils as mpl_utils

import logging
import misc.logging_utils as logging_utils
logger = logging.getLogger(__name__)

default_title = None
default_min_read_length = None
default_max_read_length = None
default_ymax = None
default_fontsize = 20

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Create a bar chart for either a single or all samples "
        "from the output of get-read-length-distribution.")
    parser.add_argument('distribution', help="The (csv) length distribution "
        "file")
    parser.add_argument('basename', help="The basename of the read counts to "
        "visualize, or ALL")
    parser.add_argument('out', help="The output (image) file")
    
    parser.add_argument('--title', help="The title of the plot", 
        default=default_title)

    parser.add_argument('--min-read-length', type=int, 
        default=default_min_read_length)
    parser.add_argument('--max-read-length', type=int, 
        default=default_max_read_length)
    parser.add_argument('--ymax', help="The maximum for the y-axis", type=int,
        default=default_ymax)
    
    parser.add_argument('--fontsize', help="The font size to use for most of "
        "the text in the plot", type=int, default=default_fontsize)

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "Reading the length distributions"
    logger.info(msg)
    distribution_df = pd.read_csv(args.distribution)

    msg = "Filtering lengths and samples based on parameters"
    logger.info(msg)

    m_sample = np.full(len(distribution_df), True)
    if args.basename != "ALL":
        m_sample = distribution_df['basename'] == args.basename
        
    m_min_read_length = np.full(len(distribution_df), True)
    if args.min_read_length is not None:
        m_min_read_length = distribution_df['length'] >= args.min_read_length

    m_max_read_length = np.full(len(distribution_df), True)
    if args.max_read_length is not None:
        m_max_read_length = distribution_df['length'] <= args.max_read_length

    m_to_view = m_min_read_length & m_max_read_length & m_sample
    distribution_df = distribution_df[m_to_view]

    # guess ymax if it is not given
    if args.ymax is None:
        msg = "Guessing ymax"
        logger.info(msg)
        args.ymax = distribution_df['count'].max()

    msg = "Creating the plots"
    logger.info(msg)

    g = sns.factorplot(
        x='length',
        y='count',
        col='basename',
        data=distribution_df,
        kind="bar",
        col_wrap=3,
        size=5,
        aspect=1.5
    )

    g.set_titles("{col_name}")
    g.set_xlabels('Length', fontsize=args.fontsize)
    g.set_ylabels('Count', fontsize=args.fontsize)

    for ax in g.axes.flat:
        # hack for punch-p data
        text = ax.title.get_text()
        
        text = text.replace(".riboseq.cell-type-hela.rep-", " (Kim, ")
        text = text.replace(".riboseq.cell-type-hela-s3.rep-", " (")
        text = text.replace(".GRCh38_85.fastq", ")")
        ax.title.set_text(text)
        
        ax.set_yscale('log')
        ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))

        ax.set_ylim((1, args.ymax))
        mpl_utils.remove_top_and_right_splines(ax)

        mpl_utils.set_title_fontsize(ax, args.fontsize)
        mpl_utils.set_ticklabels_fontsize(ax, args.fontsize)
        mpl_utils.set_ticklabel_rotation(ax, 90)
        mpl_utils.hide_first_y_tick_label(ax)
        
        
    if args.title is not None:
        g.fig.suptitle(args.title, size=args.fontsize, y=1.01)

    msg = "Writing the plot to disk"
    logger.info(msg)
    g.savefig(args.out)

if __name__ == '__main__':
    main()
