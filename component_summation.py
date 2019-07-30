#!/usr/bin/env python3

"""Plot PyBDSF componets for each source."""

import sys
import argparse
import pandas as pd
import aplpy

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '18 June 2019'


def get_ellipses(csv):
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
    pa = df['PA_2'] - 90  # PA_img_plane_2, DC_PA_2, DC_PA_img_plane_2
    flux = df['Total_flux_2']
    fraction = flux / df['Total_flux_1']

    ellipses = []
    for i, _ in enumerate(source_name):
        ellipse = (source_name[i], flux[i], fraction[i],
                   ra[i], dec[i], major[i], minor[i], pa[i])
        ellipses.append(ellipse)

    return ellipses


def get_info(csv):
    """Get the unique list of blazar names.

    Parameters
    ----------
    ellipses : list
        Tuples with the ellipse parameters.

    Returns
    -------
    list
        Names of the blazars.
    """
    df = pd.read_csv(csv)

    names = list(dict.fromkeys(df['Source name']).keys())
    ras = list(dict.fromkeys(df['RA_1']).keys())
    decs = list(dict.fromkeys(df['DEC_1']).keys())
    fluxes = list(dict.fromkeys(df['Peak_flux_1']).keys())

    return names, ras, decs, fluxes


def unpack(s):
    """Convenience function to get a list as a string without the braces.

    Parameters
    ----------
    s : list
        The list to be turned into a string.

    Returns
    -------
    string
        The list as a string.
    """
    x = " ".join(map(str, s))
    return x


def make_region_files(names, ellipses,
                      my_dir='/data5/sean/deep-fields/images/regions/'):
    """Create region files for each source.

    Parameters
    ----------
    names : list
        Blazar names.
    ellipses : list
        Tuples of the ellipse parameters.

    Returns
    -------
    list
        Names of the created region files.
    """
    header = ('# Region file format: DS9 version 4.1\nglobal color=green das' +
              'hlist=8 3 width=1 font="helvetica 10 normal roman" select=1 h' +
              'ighlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 sou' +
              'rce=1\nicrs\n')

    region_files = []
    for name in names:
        region = f'{my_dir}{name}.reg'
        region_files.append(region)
        with open(region, 'a') as the_file:
            the_file.write(header)

        for name_, flux, fraction, *ellipse in ellipses:
            if name == name_:
                ellipse = f'ellipse {unpack(ellipse)} # color=white width=2\n'
                with open(region, 'a') as the_file:
                    the_file.write(ellipse)

    return region_files


def fits_file(name, ra, my_dir='/data5/sean/deep-fields/catalogues/'):
    """Get the FITS file for a given source.

    Parameters
    ----------

    Returns
    -------
    string
        Name of the FITS file of the field the source is in.
    """

    if ra > 156 and ra < 168:
        field_file = f'{my_dir}lockman.hole.11.06.2019.img.fits'
    elif ra > 214 and ra < 222:
        field_file = f'{my_dir}bootes.11.06.2019.img.fits'
    elif ra > 237 and ra < 248:
        field_file = f'{my_dir}elias.n1.11.06.2019.img.fits'
    else:
        raise ValueError(f'{name} not in any field.')
    return field_file


def plot_ellipses(names, region_files, ras, decs, fluxes, radius=1 / 60,
                  cmap='viridis', vmin=0,
                  my_dir='/data5/sean/deep-fields/images/component-summation'):
    """Use AplPy's show_regions to plot the region files.

    Parameters
    ----------

    Returns
    -------
    """
    for name, ra, dec, flux, region_file in zip(names, ras, decs, fluxes,
                                                region_files):
        field_file = fits_file(name=name, ra=ra)
        image = aplpy.FITSFigure(field_file)  # atom gives MemoryError
        image.show_regions(region_file, set_facecolor='red')
        image.recenter(ra, dec, radius=radius)
        image.show_colorscale(cmap=cmap, vmin=vmin, vmax=flux,
                              stretch='arcsinh')
        image.save(f'{my_dir}/{name}.png')
        print(f'Done! View it: gpicview {my_dir}/{name}.png')
        sys.exit()


def main():
    """Plot PyBDSF componets for each source.

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
                        default=('/data5/sean/deep-fields/catalogues/' +
                                 'deep.fields.29.07.2019.gaul.csv'),
                        help='CSV catalogue of the blazars')

    args = parser.parse_args()

    csv = args.csv
    ellipses = get_ellipses(csv=csv)
    names, ras, decs, fluxes = get_info(csv=csv)
    region_files = make_region_files(names=names, ellipses=ellipses)
    plot_ellipses(names=names, region_files=region_files, ras=ras, decs=decs,
                  fluxes=fluxes)


if __name__ == '__main__':
    main()
