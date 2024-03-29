#!/usr/bin/env python3

"""Create a 2D cutout, and save it to a new FITS file, including the updated
cutout WCS."""

import argparse
import numpy as np
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy import wcs as PPP
from astropy.coordinates import SkyCoord
from astropy import units as u

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '18 June 2019'


def save_cutout(filename, position, size, source_name,
                my_dir='/mnt/closet/deep-fields/images/radio-fits/'):
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
    w = PPP.WCS(hdu.header, naxis=2)

    hdu.data = np.squeeze(hdu.data)  # two dimensions only

    cutout = Cutout2D(hdu.data, position=position, size=size, wcs=w)
    # hdu.header.update(cutout.wcs.to_header())  # update the fits header
    hdu.header.update(w.to_header())  # update the fits header

    cutout_filename = '/mnt/closet/deep-fields/catalogues/lockman.hole.11.06.2019.slice.fits'  # f'{my_dir}{source_name}.fits'
    hdu.writeto(cutout_filename, overwrite=True)  # save
    print(f'ds9 {cutout_filename}')


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
                        default= '/mnt/closet/deep-fields/catalogues/lockman.hole.11.06.2019.img.fits',
                        help='Bootes FITS file')

    args = parser.parse_args()
    filename = args.filename
    position = SkyCoord(161.6356597, 58.4681410, unit='deg')
    size = u.Quantity((2, 2), u.arcmin)

    save_cutout(filename=filename,
                position=position,
                size=size,
                source_name='5BZQJ1422+3223')


if __name__ == '__main__':
    main()
