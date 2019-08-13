#!/usr/bin/env python3

"""Make a few histograms and scatter plots."""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '07 August 2019'


def get_data(csv, column):
    """Get data from the CSV for a specific parameter.

    Parameters
    ----------
    csv : string
        Name of the CSV file.
    column : string
        Name of the column to be returned.

    Returns
    -------
    list
        BL Lac values for that column.
    list
        FSRQ values for that column.
    """
    df = pd.read_csv(csv)
    bllac = df[df['Type'] == 'BL Lac'][column]
    fsrq = df[df['Type'] == 'FSRQ'][column]
    return bllac, fsrq


def save_name(column):
    """Remove characters from a string to make it good for saving.

    Paramters
    ---------
    column : string

    Returns
    -------
    string
        The column name with special characters removed.
    """
    remove_these = [' (%)', ' (kpc)', ' ']
    replace_with_these = ['', '', '-']
    for remove, replace in zip(remove_these, replace_with_these):
        column = column.lower().replace(remove, replace)
    return column


def histogram(bllac, fsrq, column, x, y, n_bin):
    """Plot data for BL Lacs and FSRQs on the same histogram plot.

    Parameters
    ----------
    bllac : list
        BL Lac values.
    fsrq : list
        FSRQ values.

    Returns
    -------
    None
    """
    bins = np.linspace(0, x, n_bin)
    plt.figure(figsize=(13.92, 8.60)).patch.set_facecolor('white')

    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['mathtext.fontset'] = 'dejavuserif'
    mpl.rcParams['xtick.major.size'] = 10
    mpl.rcParams['xtick.major.width'] = 2
    mpl.rcParams['xtick.minor.size'] = 5
    mpl.rcParams['xtick.minor.width'] = 2
    mpl.rcParams['ytick.major.size'] = 10
    mpl.rcParams['ytick.major.width'] = 2
    mpl.rcParams['ytick.minor.size'] = 5
    mpl.rcParams['ytick.minor.width'] = 2
    mpl.rcParams['axes.linewidth'] = 2

    plt.hist(bllac, bins=bins, alpha=1, lw=4, histtype='step', fill=False,
             color='red', label='BL Lacs')
    plt.hist(fsrq, bins=bins, alpha=1, lw=4, histtype='step', fill=False,
             color='blue', label='FSRQs', linestyle='--')

    plt.xlabel(column, fontsize=30, color='black')
    plt.ylabel('$N$', fontsize=30, color='black')
    plt.xlim(0, x)
    plt.ylim(0, y)
    plt.xticks(np.linspace(0, x, n_bin), fontsize=30, color='black')
    plt.yticks(np.linspace(0, y, y + 1), fontsize=30, color='black')

    # ax = plt.gca()
    # yticks = ax.yaxis.get_major_ticks()
    # yticks[0].label1.set_visible(False)

    legend = plt.legend(bbox_to_anchor=(0, 1.0, 1, 0), loc='lower left',
                        mode='expand', numpoints=1, fontsize=30, ncol=2,
                        frameon=False)
    plt.setp(legend.get_texts(), color='black')

    column = save_name(column=column)
    plt.savefig(f'/mnt/closet/deep-fields/results/{column}.png',
                bbox_inches='tight', format='png')
    # plt.show()
    plt.close()


def scatter(bllac_x, bllac_y, fsrq_x, fsrq_y, column_x, column_y, xlim, ylim):
    """Scatter plot for two quantities for BL Lacs and FSRQs.

    Parameters
    ----------
    bllac_x : list
        BL Lac data for the X axis.
    bllac_y : list
        BL Lac data for the Y axis.
    fsrq_x : list
        FSRQ data for the X axis.
    fsrq_y : list
        FSRQ data for the Y axis.
    Returns
    -------
    None
    """

    plt.figure(figsize=(13.92, 8.60)).patch.set_facecolor('white')

    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['mathtext.fontset'] = 'dejavuserif'
    mpl.rcParams['xtick.major.size'] = 10
    mpl.rcParams['xtick.major.width'] = 2
    mpl.rcParams['xtick.minor.size'] = 5
    mpl.rcParams['xtick.minor.width'] = 2
    mpl.rcParams['ytick.major.size'] = 10
    mpl.rcParams['ytick.major.width'] = 2
    mpl.rcParams['ytick.minor.size'] = 5
    mpl.rcParams['ytick.minor.width'] = 2
    mpl.rcParams['axes.linewidth'] = 2

    plt.plot(bllac_x, bllac_y, marker='o', color='red', ls='None',
             label='BL Lacs', ms=10)
    plt.plot(fsrq_x, fsrq_y, marker='o', color='blue', ls='None',
             label='FSRQs', ms=10)

    plt.xlabel(column_x, fontsize=30, color='black')
    plt.ylabel(column_y, fontsize=30, color='black')
    plt.xlim(xlim[0] - 5, xlim[1])
    plt.ylim(ylim[0], ylim[1])
    plt.xticks(np.linspace(xlim[0], xlim[1], xlim[2]), fontsize=30,
               color='black')
    plt.yticks(np.linspace(ylim[0], ylim[1], ylim[2]), fontsize=30,
               color='black')

    legend = plt.legend(bbox_to_anchor=(0, 1.0, 1, 0), loc='lower left',
                        mode='expand', numpoints=1, fontsize=30, ncol=2,
                        frameon=False)
    plt.setp(legend.get_texts(), color='black')

    column_x = save_name(column=column_x)
    column_y = save_name(column=column_y)
    plt.savefig(f'/mnt/closet/deep-fields/results/{column_x}-v-{column_y}.png',
                bbox_inches='tight', format='png')
    # plt.show()
    plt.close()


def main():
    """Plot quantities.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    formatter_class = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=formatter_class)
    parser.add_argument('-c',
                        '--csv',
                        required=False,
                        type=str,
                        default='/home/sean/Downloads/ldf1-results - results' +
                        '.csv',
                        help='CSV catalogue of the blazars')

    args = parser.parse_args()
    csv = args.csv
    column1 = 'Diffuse emission (%)'
    bllac_diffuse, fsrq_diffuse = get_data(csv=csv, column=column1)
    histogram(bllac=bllac_diffuse, fsrq=fsrq_diffuse, column=column1, x=100,
              y=7, n_bin=11)

    column2 = 'Extent (kpc)'
    bllac_extent, fsrq_extent = get_data(csv=csv, column=column2)
    scatter(bllac_x=bllac_diffuse, bllac_y=bllac_extent, fsrq_x=fsrq_diffuse,
            fsrq_y=fsrq_extent, column_x=column1, column_y=column2,
            xlim=[0, 100, 11], ylim=[50, 400, 8])


if __name__ == '__main__':
    main()
