#!/usr/bin/env python3

"""Create a 2D cutout, and save it to a new FITS file, including the updated
cutout WCS."""

import argparse
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '18 June 2019'


def save_cutout(filename, position, size, source_name,
                my_dir='/data5/sean/deep-fields/images/radio-fits/'):
    """Make postage stamp and save it.

    Parameters
    ----------
    filename : string
        The FITS filename that is to be cropped.
    position : tuple
        The X and Y position that the cropped image is to be centred on.
    size : tuple
        The size of the final cropped image.

    Returns
    -------
    None
    """
    hdu = fits.open(filename)[0]  # load the image and the wcs
    wcs = WCS(hdu.header)

    # make the cutout, including the wcs
    cutout = Cutout2D(hdu.data, position=position, size=size, wcs=wcs)

    # Put the cutout image in the fits hdu
    hdu.data = cutout.data

    # update the fits header with the cutout wcs
    hdu.header.update(cutout.wcs.to_header())

    # Write the cutout to a new FITS file
    cutout_filename = f'{my_dir}{source_name}.fits'
    hdu.writeto(cutout_filename, overwrite=True)


def main():
    """Create a 2D cutout, and save it to a new FITS file, including the
    updated cutout WCS.
    """
    formatter_class = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=formatter_class)

    parser.add_argument('-f',
                        '--filename',
                        required=False,
                        type=str,
                        default='/data5/sean/deep-fields/catalogues/bootes.11.06.2019.img.fits',
                        help='Bootes FITS file')

    args = parser.parse_args()
    filename = args.filename
    # position = SkyCoord(215.62658, 32.38622, unit='deg')
    # size = u.Quantity((2, 2), u.arcmin)

    save_cutout(filename=filename,
                position=[215.62658*u.deg, 32.38622*u.deg, 0, 0], size=[2,2,:,:], source_name='5BZQJ1422+3223')


if __name__ == '__main__':
    main()
