#!/usr/bin/env python3

"""Measure the extent of sources."""

import sys
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from ds9norm import DS9Normalize
from scipy import ndimage

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '06 August 2019'


def abline(slope, intercept):
    """Plot a line from slope and intercept."""
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, '--')


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

    if name == 'ILTJ125524.94+544942.7':
        data[90:, :] = np.nan
        data[:35, :] = np.nan

    if name == 'ILTJ125544.98+553748.6':
        data[:35, :] = np.nan
        data[63:71, 32:39] = np.nan

    if name == 'ILTJ125624.63+552823.7':
        data[:, 88:] = np.nan

    if name == 'ILTJ125755.66+541155.1':
        d[:, :37] = np.nan
        d[:, 83:] = np.nan

    return data


my_dir = '/home/sean/Downloads/'  # /data5/sean
df = pd.read_csv(f'{my_dir}fr-i-fr-ii-mingo-ldr1.csv')  # load data
df = df[df['Mosaic_ID'] == 'P196+55']  # just doing P196+55 for now
fields = [f'{my_dir}{m}-mosaic.fits' for m in df['Mosaic_ID']]
thresholds = 5 * df['Isl_rms'] / 1000  # converting to jansky
plt.figure(figsize=(8, 8))
i = 0

for source_name, ra, dec, field, threshold in zip(df['Source_Name_1'],
                                                  df['RA_1'],
                                                  df['DEC_1'], fields,
                                                  thresholds):
    # if source_name != '5BZQJ1604+5714':  # for testing one source
    #     continue
    hdu = fits.open(field)[0]
    wcs = WCS(hdu.header, naxis=2)
    sky_position = SkyCoord(ra, dec, unit='deg')
    cutout = Cutout2D(np.squeeze(hdu.data), sky_position,
                      size=[3, 3] * u.arcmin, wcs=wcs)

    d = cutout.data
    d[d < threshold] = np.nan
    d_plot = d[:]

    # set inside elements to nan to speed up the calculations
    for x in range(0, d.shape[0]):
        for y in range(0, d.shape[0]):
            try:
                if (np.isfinite(d[x - 1, y]) and np.isfinite(d[x + 1, y]) and
                   np.isfinite(d[x, y - 1]) and np.isfinite(d[x, y + 1])):
                    d[x, y] = -3.141592  # can use any number
            except IndexError:
                print(f'IndexError encountered for {source_name}.')
            except:
                raise

    d[d == -3.141592] = np.nan
    d = manual_mask(name=source_name, data=d)
    d = d / d
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

    slope = (max_y2 - max_y1) / (max_x2 - max_x1)

    # for the good cells, find the distance from it to the point on the line
    # we assume we are looking down the line of longest length and the source
    # is flat in the plane of the sky
    # https://stackoverflow.com/a/52756183/6386612

    p1 = np.array([max_x1, max_y1])
    p2 = np.array([max_x2, max_y2])
    widths = []
    width_x, width_y = [], []

    for (x, y) in good_cells:
        p3 = np.array([x, y])
        widths.append(np.cross(p2 - p1, p3 - p1) / np.linalg.norm(p2 - p1))
        width_x.append(x)
        width_y.append(y)

    widths = np.array(widths)
    my_max_width = np.max(widths)
    my_min_width = np.min(widths)
    width_x_max = width_x[widths.argmax()]
    width_y_max = width_y[widths.argmax()]
    width_x_min = width_x[widths.argmin()]
    width_y_min = width_y[widths.argmin()]
    asec_width = (my_max_width + np.abs(my_min_width)) * 1.5  # 1.5" per pixel

    fig = plt.imshow(d, vmin=0, vmax=np.nanmax(d), origin='lower',
                     cmap='Greys')  # ,norm=DS9Normalize(stretch='arcsinh'))
    # fig = plt.imshow(d, vmin=0, vmax=np.nanmax(d), origin='lower',
    #                  norm=DS9Normalize(stretch='arcsinh'))
    plt.colorbar()
    fig.axes.get_xaxis().set_visible(False)
    fig.axes.get_yaxis().set_visible(False)
    circle = plt.Circle(centre, my_max / 2, color='r', fill=False, alpha=0.5,
                        lw=2)
    fig = plt.gcf()
    ax = fig.gca()
    # ax.add_artist(circle)

    plt.plot([max_y1, max_y2], [max_x1, max_x2], color='red', alpha=0.5, lw=2)
    plt.plot([width_y_max], [width_x_max], marker='o', markersize=10,
             color='red', alpha=0.5)
    plt.plot([width_y_min], [width_x_min], marker='o', markersize=10,
             color='red', alpha=0.5)

    theta = math.radians(90)
    cx = (max_x1 + max_x2) / 2
    cy = (max_y1 + max_y2) / 2
    x1_ = (((max_x1 - cx) * math.cos(theta) + (max_y1 - cy) * math.sin(theta))
           + cx)
    x2_ = (((max_x2 - cx) * math.cos(theta) + (max_y2 - cy) * math.sin(theta))
           + cx)
    y1_ = ((-(max_x1 - cx) * math.sin(theta) + (max_y1 - cy) * math.cos(theta))
           + cy)
    y2_ = ((-(max_x2 - cx) * math.sin(theta) + (max_y2 - cy) * math.cos(theta))
           + cy)

    # plt.plot([y1_, y2_], [x1_, x2_], color='green', alpha=0.5, lw=2)
    # plt.plot([(max_y1 + max_y2) / 2], [(max_x1 + max_x2) / 2], marker='o',
    #          color='red', alpha=0.5)

    plt.title(f'{source_name}\n5\u03C3 = {threshold * 1000:.3f} mJy; d = ' +
              f'width = {asec_width:.1f}"')
    plt.show()
    # plt.savefig(f'{my_dir}../images/extention/{source_name}.png')
    plt.clf()
    print(f'{source_name} {asec_max} {asec_width} {threshold * 1000}')

    i += 1
    if i > 8:
        sys.exit()
