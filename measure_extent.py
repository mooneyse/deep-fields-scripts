#!/usr/bin/env python3

"""Measure the extent of sources."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from ds9norm import DS9Normalize

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '06 August 2019'


def manual_mask(name, data):
    """Some ad hoc adjustments to a few sources, made after visual inspection.

    Parameters
    ----------
    name : string
        Source name.
    data : array
        pixel values.

    Returns
    -------
    array
        Values with the manual masking done.
    """
    if name == '5BZQJ1422+3223':
        d[:, 70:] = np.nan
        d[20:23, 51:54] = np.nan
        d[42:45, 58:60] = np.nan

    if name == '5BZBJ1426+3404':
        data[:29, :] = np.nan
        data[63:, :] = np.nan

    if name == '5BZQJ1429+3529':
        data[:, :20] = np.nan
        data[:26, :] = np.nan

    if name == '5BZQJ1435+3353':
        data[:, :15] = np.nan
        data[:, 47:] = np.nan

    if name == '5BZQJ1437+3618':
        data[:, :23] = np.nan

    if name == '5BZQJ1437+3519':
        data[:, 65:] = np.nan

    if name == '5BZBJ1558+5625':
        data[62:, :] = np.nan

    if name == '5BZBJ1605+5421':
        data[:, 60:] = np.nan
        data[:, :29] = np.nan
        data[66:, :] = np.nan
        data[:31, :] = np.nan

    if name == '5BZQJ1606+5405':
        data[:, :20] = np.nan
        data[:, 62:] = np.nan

    if name == '5BZQJ1608+5613':
        data[:30, :] = np.nan
        data[:, 47:] = np.nan
        data[46:54, 43:47] = np.nan

    if name == '5BZQJ1619+5256':
        data[77:, :] = np.nan
        data[:32, :] = np.nan
        data[43:47, 46:54] = np.nan

    if name == '5BZBJ1037+5711':
        data[:25, :] = np.nan
        data[54:59, 37:44] = np.nan

    if name == '5BZQJ1604+5714':
        data[38:49, 22:30] = np.nan

    return data


my_dir = '/mnt/closet/deep-fields/catalogues/'  # /data5/sean
df = pd.read_csv(f'{my_dir}deep.fields.11.06.2019.cat.csv')  # load data
fields = [f'{my_dir}{f}.11.06.2019.img.fits'.lower().replace(' ', '.')
          for f in df['Field']]
thresholds = 5 * df['Isl_rms']
plt.figure(figsize=(8, 8))

for source_name, ra, dec, field, threshold in zip(df['Source name'], df['RA'],
                                                  df['DEC'], fields,
                                                  thresholds):
    # if source_name != '5BZQJ1604+5714':  # for testing one source
    #     continue
    hdu = fits.open(field)[0]
    wcs = WCS(hdu.header, naxis=2)
    sky_position = SkyCoord(ra, dec, unit='deg')
    cutout = Cutout2D(np.squeeze(hdu.data), sky_position,
                      size=[2, 2] * u.arcmin, wcs=wcs)

    d = cutout.data
    d[d < threshold] = np.nan
    d = manual_mask(name=source_name, data=d)
    rows, cols = d.shape
    good_cells = []

    for r in range(rows):
        for c in range(cols):
            if not np.isnan(d[r, c]):
                good_cells.append([r, c])

    # find distance between good_cell and all other good_cells
    distances, x1s, x2s, y1s, y2s = [], [], [], [], []
    for (x1, y1) in good_cells:
        for (x2, y2) in good_cells:
            distances.append(np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2))
            x1s.append(x1)
            x2s.append(x2)
            y1s.append(y1)
            y2s.append(y2)

    distances = np.array(distances)
    my_max = np.max(distances)
    max_x1 = x1s[distances.argmax()]
    max_x2 = x2s[distances.argmax()]
    max_y1 = y1s[distances.argmax()]
    max_y2 = y2s[distances.argmax()]
    asec_max = my_max * 1.5  # 1.5" per pixel
    centre = ((max_y1 + max_y2) / 2, (max_x1 + max_x2) / 2)

    fig = plt.imshow(d, vmin=0, vmax=np.nanmax(d), origin='lower',
                     norm=DS9Normalize(stretch='arcsinh'))
    plt.colorbar()
    fig.axes.get_xaxis().set_visible(False)
    fig.axes.get_yaxis().set_visible(False)
    circle = plt.Circle(centre, my_max / 2, color='r', fill=False, alpha=0.5,
                        lw=2)
    fig = plt.gcf()
    ax = fig.gca()
    ax.add_artist(circle)

    plt.plot([max_y1, max_y2], [max_x1, max_x2], color='red', alpha=0.5, lw=2)
    plt.plot([(max_y1 + max_y2) / 2], [(max_x1 + max_x2) / 2], marker='o',
             color='red', alpha=0.5)
    plt.title(f'{source_name}\n5\u03C3 = {threshold * 1000:.3f} mJy; d = ' +
              f'{asec_max:.1f}"')
    # plt.show()
    plt.savefig(f'{my_dir}../images/extention/{source_name}.png')
    plt.clf()
    print(f'{source_name} {asec_max} {threshold * 1000}')
