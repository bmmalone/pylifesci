#! /usr/bin/env python3
import matplotlib
matplotlib.use('agg')

import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import misc.mpl_utils as mpl_utils

import logging
import misc.logging_utils as logging_utils
logger = logging.getLogger(__name__)

default_y_tick_frequency = 1e7
default_title = None

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script creates a stacked bar chart of the uniquely- "
        "and multi-mapping read counts found by count-aligned-reads")

    parser.add_argument('read_counts', help="The counts csv file from "
        "count-aligned-reads")
    parser.add_argument('out', help="The output image file")

    parser.add_argument('--y-tick-frequency', help="The frequency of the ticks "
        "on the y-axis", type=int, default=default_y_tick_frequency)

    parser.add_argument('-t', '--title', help="The title for the plot",
        default=default_title)

    
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "Reading count"
    logger.info(msg)
    read_counts = pd.read_csv(args.read_counts)

    msg = "Extracting bar heights"
    logger.info(msg)

    fields = ['uniquely_aligned_reads', 'multimapping_reads']
    pretty_fields = ['Unique', 'Multimapping']
    read_counts_np = read_counts[fields].values

    msg = "Determing y-ticks"
    logger.info(msg)

    ytm = np.max(read_counts['aligned_reads']) / args.y_tick_frequency
    y_tick_max = round(ytm) * args.y_tick_frequency
    y_ticks = np.arange(0, y_tick_max, args.y_tick_frequency)
    y_tick_labels = ["{:.1e}".format(y) for y in y_ticks]

    msg = "Plotting the bars"
    logger.info(msg)

    x_tick_labels = read_counts['dataset']

    fig, ax = plt.subplots()

    mpl_utils.create_stacked_bar_plot(ax, read_counts_np, 
                                      stack_labels=pretty_fields,
                                      colors=plt.cm.Blues, edge_colors='k', 
                                      gap=0.1, end_gaps=True, 
                                      y_ticks=y_ticks, y_tick_labels=y_tick_labels,
                                      x_tick_labels=x_tick_labels
                                     )

    if args.title is not None:
        ax.set_title(args.title)

    msg = "Saving plot to disk"
    logger.info(msg)

    fig.savefig(args.out, bbox_inches='tight')


if __name__ == '__main__':
    main()
