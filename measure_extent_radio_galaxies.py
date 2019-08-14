#!/usr/bin/env python3

"""Measure the extent of sources."""

import matplotlib as mpl
mpl.use('Agg')
import os.path
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

    if name == 'ILTJ125805.10+541704.4':
        d[97:, :] = np.nan
        d[:, 88:] = np.nan
        d[:, :20] = np.nan

    if name == 'ILTJ125840.41+534935.5':
        d[80:, :] = np.nan

    if name == 'ILTJ125852.58+555527.3':
        d[15:42, :18] = np.nan

    if name == 'ILTJ125959.20+560735.0':
        d[100:, :] = np.nan
        d[:, 80:] = np.nan

    if name == 'ILTJ130005.50+551315.1':
        d[:, :14] = np.nan

    if name == 'ILTJ130018.09+545048.2':
        d[:, 80:] = np.nan

    if name == 'ILTJ130057.57+542654.6':
        d[29:37, 50:60] = np.nan

    if name == 'ILTJ130102.49+561003.6':
        d[100:, :] = np.nan
        d[:30, :] = np.nan
        d[61:68, 80:87] = np.nan

    if name == 'ILTJ130112.93+550410.5':
        d[90:, :] = np.nan

    if name == 'ILTJ130138.23+551918.4':
        d[:, 83:] = np.nan

    if name == 'ILTJ130148.82+544723.9':
        d[114:, :] = np.nan

    if name == 'ILTJ130152.93+555557.4':
        d[:, :30] = np.nan
        d[:, 96:] = np.nan

    if name == 'ILTJ130159.04+534548.8':
        d[:34, :] = np.nan
        d[103:, :] = np.nan

    if name == 'ILTJ130528.93+561918.9':
        d[:, 90:] = np.nan

    if name == 'ILTJ130605.63+555127.6':
        d[:, 96:] = np.nan

    if name == 'ILTJ130634.72+553657.7':
        d[107:, :] = np.nan

    if name == 'ILTJ130804.04+550835.4':
        d[:, :46] = np.nan
        d[:, 70:] = np.nan
        d[:13, :] = np.nan

    if name == 'ILTJ130821.80+562744.1':
        d[:, :23] = np.nan
        d[:, 100:] = np.nan

    if name == 'ILTJ130857.14+542915.9':
        d[:, 98:] = np.nan
        d[:, :23] = np.nan

    if name == 'ILTJ130907.51+545520.4':
        d[:, 96:] = np.nan

    if name == 'ILTJ130926.29+534820.0':
        d[54:57, 76:79] = np.nan
        d[44:53, 69:74] = np.nan

    if name == 'ILTJ131115.53+534356.8':
        d[101:, :] = np.nan

    if name == 'ILTJ131236.04+535524.7':
        d[:, 83:] = np.nan

    if name == 'ILTJ125757.12+554224.1':
        d[:, :40] = np.nan
        d[:, 71:] = np.nan

    if name == 'ILTJ125840.51+534935.5':
        d[79:, :] = np.nan
        # d[:45, :] = np.nan
        d[:, 70:] = np.nan

    if name == 'ILTJ125852.58+555527.3':
        d[:, :21] = np.nan
        d[71:, :] = np.nan

    if name == 'ILTJ130005.50+551315.1':
        d[:, :32] = np.nan

    if name == 'ILTJ130112.93+550410.5':
        d[:33, :] = np.nan

    if name == 'ILTJ130148.82+544723.9':
        d[:50, :] = np.nan

    if name == 'ILTJ130331.44+534400.9':
        d[:, :43] = np.nan
        d[67:85, 63:88] = np.nan
        d[71:, :] = np.nan
        d[:, 72:] = np.nan
        d[65:85, 66:88] = np.nan

    if name == 'ILTJ130441.11+551952.7':
        d[:, :39] = np.nan
        d[:, 73:] = np.nan
        d[:37, :] = np.nan

    if name == 'ILTJ130531.20+535432.3':
        d[:, :39] = np.nan
        d[:, 97:] = np.nan
        d[:30, :] = np.nan
        d[:30, :] = np.nan
        d[45:53, 55:64] = np.nan

    if name == 'ILTJ130605.63+555127.6':
        d[:, :41] = np.nan

    if name == 'ILTJ130634.72+553657.7':
        d[73:85, 25:42] = np.nan
        d[30:41, 43:56] = np.nan

    if name == 'ILTJ130638.74+541448.8':
        d[:, :22] = np.nan
        d[:, 91:] = np.nan
        d[:45, :] = np.nan

    if name == 'ILTJ130749.04+545120.0':
        d[:, :11] = np.nan
        d[91:, :] = np.nan
        d[:42, :] = np.nan
        d[42:51, 54:63] = np.nan

    if name == 'ILTJ130804.04+550836.4':
        d[:17, :] = np.nan

    if name == 'ILTJ130857.14+542915.9':
        d[58:68, 35:45] = np.nan
        d[66:77, 35:52] = np.nan
        d[:, 71:] = np.nan
        d[40:49, 60:76] = np.nan

    if name == 'ILTJ130926.29+534820.0':
        d[65:77, 42:55] = np.nan

    if name == 'ILTJ130331.30+540024.1':
        d[174:, :] = np.nan
        d[:, :50] = np.nan

    return data


my_dir = '/data5/sean/deep-fields/catalogues/'  # mnt/closet
df = pd.read_csv(f'{my_dir}radio.galaxies.ldr1.csv')  # load data
# df = df[df['Mosaic_ID'] == 'P196+55']  # just doing P196+55 for now
fields = [f'/data1/lotss-data/{m}-mosaic.fits' for m in df['Mosaic_ID']]
thresholds = 5 * df['Isl_rms'] / 1000  # converting to jansky
plt.figure(figsize=(8, 8))

results_csv = f'{my_dir}../results/extention-radio-galaxies.csv'
result_header = ('name,type,LM flux,LM size,length,width,threshold,z,WAT,' +
                 'NAT,D-D\n')

if not os.path.exists(results_csv):
    with open(results_csv, 'a') as f:
        f.write(result_header)

i = 0

for (source_name, ra, dec, field, threshold, fri, frii, redshift, nat, wat,
     dd, lm_size, lm_flux) in zip(df['Source_Name_1'], df['RA_1'],
                                  df['DEC_1'], fields, thresholds,
                                  df['FR1'], df['FR2'], df['z_best'],
                                  df['NAT'], df['WAT'], df['D-D'],
                                  df['LM_dec_size'], df['LM_Flux']):

    save = f'{my_dir}../images/extention-radio-galaxies/{source_name}.png'
    if os.path.exists(save):
        print(f'{save} already exists.')
        continue

    source_type = 'FR-I' if fri else 'FR-II'  # already removed small sources

    # if source_name != 'ILTJ130804.04+550835.4':  # for testing one source
    #     continue

    try:
        hdu = fits.open(field)[0]
    except FileNotFoundError:  # fits might have different name
        print(f'{field} does not exist.')
        continue

    wcs = WCS(hdu.header, naxis=2)
    sky_position = SkyCoord(ra, dec, unit='deg')
    s = math.ceil((lm_size * 2) / 60)
    size = [3, 3] if field[-19:-5] == 'P196+55-mosaic' else [s, s]

    # some sources need a bigger cutout
    if source_name == 'ILTJ130140.02+540825.9':
        size = [4, 4]
    elif source_name == 'ILTJ130331.30+540024.1':
        size = [6, 6]

    cutout = Cutout2D(np.squeeze(hdu.data), sky_position,
                      size=size * u.arcmin, wcs=wcs)

    d = cutout.data
    d[d < threshold] = np.nan

    # set inside elements to nan to speed up the calculations
    # for x in range(0, d.shape[0]):
    #     for y in range(0, d.shape[0]):
    #         try:
    #             if (np.isfinite(d[x - 1, y]) and np.isfinite(d[x + 1, y]) and
    #                np.isfinite(d[x, y - 1]) and np.isfinite(d[x, y + 1])):
    #                 d[x, y] = -3.141592  # can use any number
    #         except IndexError:
    #             print(f'IndexError encountered for {source_name}.')
    #         except Exception:  # raise other errors as usual
    #             raise
    #
    # d[d == -3.141592] = np.nan  # have to do it this way to avoid a pattern
    d = manual_mask(name=source_name, data=d)
    # d = d / d  # scale to unity
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
    try:
        my_max = np.max(distances)
    except ValueError:
        print(f'Failed for {source_name}.')
        continue

    max_x1 = x1s[distances.argmax()]
    max_x2 = x2s[distances.argmax()]
    max_y1 = y1s[distances.argmax()]
    max_y2 = y2s[distances.argmax()]
    asec_max = my_max * 1.5  # 1.5" per pixel
    centre = ((max_y1 + max_y2) / 2, (max_x1 + max_x2) / 2)

    # for good cells find the distance from it to the point on the line, we
    # assume we look down the longest line and the source is flat in the plane
    # of the sky https://stackoverflow.com/a/52756183/6386612

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

    # fig = plt.imshow(d, vmin=0, vmax=np.nanmax(d), origin='lower',
    #                  cmap='Greys')
    fig = plt.imshow(d, vmin=0, vmax=np.nanmax(d), origin='lower',
                     norm=DS9Normalize(stretch='arcsinh'))

    plt.colorbar()
    # fig.axes.get_xaxis().set_visible(False)
    # fig.axes.get_yaxis().set_visible(False)

    plt.plot([max_y1, max_y2], [max_x1, max_x2], color='red', alpha=0.5, lw=2)
    plt.plot([width_y_max], [width_x_max], marker='o', markersize=10,
             color='red', alpha=0.5)
    plt.plot([width_y_min], [width_x_min], marker='o', markersize=10,
             color='red', alpha=0.5)

    plt.minorticks_on()
    # plt.gca().xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(10))
    # plt.gca().yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(10))
    plt.grid(which="both", linewidth=0.72, color="k")
    plt.tick_params(which="minor", length=0)

    my_string = ''
    if nat:
        my_string += 'NAT '
    elif wat:
        my_string += 'WAT '
    if dd:
        my_string += 'D-D '

    plt.title(f'{source_name}\n{source_type}; {my_string}5\u03C3 = ' +
              f'{threshold * 1000:.3f} mJy; width = {asec_width:.1f}"')
    # plt.show()
    plt.savefig(save)
    plt.clf()
    results = (f'{source_name},{source_type},{lm_flux},{lm_size},{asec_max},' +
               f'{asec_width},{threshold * 1000},{redshift},{nat},{wat},{dd}' +
               '\n')

    with open(results_csv, 'a') as f:
        f.write(results)

    i += 1
    if i > 100:
        import sys
        sys.exit()
    # build clean sample on visual inspection
