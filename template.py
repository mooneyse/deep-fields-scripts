#!/usr/bin/env python3

"""asdf"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '18 June 2019'


def asdf(sean=''):
    """Short summary.

    Parameters
    ----------
    p : type
        Description of parameter `p`.

    Returns
    -------
    type
        Description of returned object.

    """
    np.sqrt(2)
    pd.read_csv(sean)
    plt.show(sean)
    return sean


def main():
    """Short summary.

    Parameters
    ----------
    asdf : string

    Returns
    -------
    type
        Description of returned object.
    """

    formatter_class = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=formatter_class)
    parser.add_argument('-c',
                        '--csv',
                        required=False,
                        type=str,
                        default='csv.csv',
                        help='CSV catalogue of the blazars')

    args = parser.parse_args()
    csv = args.csv
    asdf(sean=csv)


if __name__ == '__main__':
    main()
