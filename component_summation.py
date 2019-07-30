#!/usr/bin/env python3

"""Plot PyBDSF componets for each source."""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import aplpy

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '18 June 2019'


def get_ellipse(csv=''):
    """Get the ellipse parameters from the CSV.

    Parameters
    ----------
    csv : string
        Filepath of the CSV containing the Gaussians fitted by PyBDSF.

    Returns
    -------
    tuple
        Associated source name, right ascension, declination, semi-major axis,
        semi-minor axis, position angle, total flux density (Jy) and fractional
        flux density (Jy) for each row in the CSV.
    """
    df = pd.read_csv(csv)

    source_name = df['Source name']
    ra = df['RA_2']
    dec = df['DEC_2']
    major = df['Maj_2']  # Maj_img_plane_2, DC_Maj_2, DC_Maj_img_plane_2
    minor = df['Min_2']  # Min_img_plane_2, DC_Min_2, DC_Min_img_plane_2
    pa = df['PA_2']  # PA_img_plane_2, DC_PA_2, DC_PA_img_plane_2
    flux = df['Total_flux_2']
    fraction = flux / df['Total_flux_1']

    ellipses = []
    for i, _ in enumerate(source_name):
        ellipse = (source_name[i], flux[i], fraction[i],
                   ra[i], dec[i], major[i], minor[i], pa[i])
        ellipses.append(ellipse)

    return ellipses


def plot_ellipses():
    """Use AplPy's show_regions to plot the region files.

    """


def main():
    """Plot PyBDSF componets for each source.

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
                        default=('/mnt/closet/deep-fields/catalogues/' +
                                 'deep.fields.29.07.2019.gaul.csv'),
                        help='CSV catalogue of the blazars')

    args = parser.parse_args()

    csv = args.csv
    ellipses = get_ellipse(csv=csv)

    header = ('# Region file format: DS9 version 4.1\nglobal color=green das' +
              'hlist=8 3 width=1 font="helvetica 10 normal roman" select=1 h' +
              'ighlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 sou' +
              'rce=1\nicrs\n')

    names = set()  # unique set of blazars
    for name, *_ in ellipses:
        names.add(name)

    def unpack(s):
        return " ".join(map(str, s))  # map(), just for kicks

    for name in names:
        region = f'/mnt/closet/deep-fields/images/regions/{name}.reg'
        with open(region, 'a') as the_file:
            the_file.write(header)

        for name_, flux, fraction, *ellipse in ellipses:
            if name == name_:
                ellipse = f'ellipse {unpack(ellipse)} # color=white width=2\n'
                with open(region, 'a') as the_file:
                    the_file.write(ellipse)

                # fits = f'/mnt/closet/deep-fields/images/sdss-fits/sdss-{name}.fits'
                # f = aplpy.FITSFigure(fits)
                # f.show_grayscale()
                # plt.show()
                # plot_ellipse(ellipse)

        # import sys
        # sys.exit()


if __name__ == '__main__':
    main()
