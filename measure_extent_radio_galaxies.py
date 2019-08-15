#!/usr/bin/env python3

"""Measure the extent of sources."""

import matplotlib as mpl
mpl.use('Agg')
import os.path
from math import sqrt, exp, sin
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


def get_dl_and_kpc_per_asec(z, H0=70, WM=0.26, WV=0.74):
    """Ned Wright's cosmology calculator.

    Returns
    -------
    float
        Luminosity distance in metres.
    float
        Kiloparsecs per arcsecond.
    """
    WR = 0.  # Omega(radiation)
    WK = 0.  # Omega curvaturve = 1-Omega(total)
    c = 299792.458  # velocity of light in km/sec
    DTT = 0.5  # time from z to now in units of 1/H0
    age = 0.5  # age of Universe in units of 1/H0
    zage = 0.1  # age of Universe at redshift z in units of 1/H0
    DCMR = 0.0  # comoving radial distance in units of c/H0
    DA = 0.0  # angular size distance
    DA_Mpc = 0.0
    kpc_DA = 0.0
    DL = 0.0  # luminosity distance
    DL_Mpc = 0.0
    a = 1.0  # 1/(1+z), the scale factor of the Universe
    az = 0.5  # 1/(1+z(object))
    h = H0 / 100.
    WR = 4.165E-5 / (h * h)  # includes 3 massless neutrino species,
    WK = 1 - WM - WR - WV
    az = 1.0 / (1 + 1.0 * z)
    age = 0.
    n = 1000  # number of points in integrals
    for i in range(n):
        a = az * (i + 0.5) / n
        adot = sqrt(WK + (WM / a) + (WR / (a * a)) + (WV * a * a))
        age = age + 1. / adot
    zage = az * age / n
    DTT = 0.0
    DCMR = 0.0
    for i in range(n):
        a = az + (1 - az) * (i + 0.5) / n
        adot = sqrt(WK + (WM / a) + (WR / (a * a)) + (WV * a * a))
        DTT = DTT + 1. / adot
        DCMR = DCMR + 1. / (a * adot)
    DTT = (1. - az) * DTT / n
    DCMR = (1. - az) * DCMR / n
    age = DTT + zage
    ratio = 1.00
    x = sqrt(abs(WK)) * DCMR
    if x > 0.1:
        if WK > 0:
            ratio = 0.5 * (exp(x) - exp(-x)) / x
        else:
            ratio = sin(x) / x
    else:
        y = x * x
        if WK < 0:
            y = -y
        ratio = 1. + y / 6. + y * y / 120.
    DCMT = ratio * DCMR
    DA = az * DCMT
    DA_Mpc = (c / H0) * DA
    kpc_DA = DA_Mpc / 206.264806
    DL = DA / (az * az)
    DL_Mpc = (c / H0) * DL
    ratio = 1.00
    x = sqrt(abs(WK)) * DCMR
    if x > 0.1:
        if WK > 0:
            ratio = (0.125 * (exp(2. * x) - exp(-2. * x)) - x / 2.) / (
                    x * x * x / 3.)
        else:
            ratio = (x / 2. - sin(2. * x) / 4.) / (x * x * x / 3.)
    else:
        y = x * x
        if WK < 0:
            y = -y
        ratio = 1. + y / 5. + (2. / 105.) * y * y
    return DL_Mpc * 3.086e22, kpc_DA


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
        data[:, :37] = np.nan
        data[:, 83:] = np.nan

    if name == 'ILTJ125805.10+541704.4':
        data[97:, :] = np.nan
        data[:, 88:] = np.nan
        data[:, :20] = np.nan

    if name == 'ILTJ125840.41+534935.5':
        data[80:, :] = np.nan

    if name == 'ILTJ125852.58+555527.3':
        data[15:42, :18] = np.nan

    if name == 'ILTJ125959.20+560735.0':
        data[100:, :] = np.nan
        data[:, 80:] = np.nan

    if name == 'ILTJ130005.50+551315.1':
        data[:, :14] = np.nan

    if name == 'ILTJ130018.09+545048.2':
        data[:, 80:] = np.nan

    if name == 'ILTJ130057.57+542654.6':
        data[29:37, 50:60] = np.nan

    if name == 'ILTJ130102.49+561003.6':
        data[100:, :] = np.nan
        data[:30, :] = np.nan
        data[61:68, 80:87] = np.nan

    if name == 'ILTJ130112.93+550410.5':
        data[90:, :] = np.nan

    if name == 'ILTJ130138.23+551918.4':
        data[:, 83:] = np.nan

    if name == 'ILTJ130148.82+544723.9':
        data[114:, :] = np.nan

    if name == 'ILTJ130152.93+555557.4':
        data[:, :30] = np.nan
        data[:, 96:] = np.nan

    if name == 'ILTJ130159.04+534548.8':
        data[:34, :] = np.nan
        data[103:, :] = np.nan

    if name == 'ILTJ130528.93+561918.9':
        data[:, 90:] = np.nan

    if name == 'ILTJ130605.63+555127.6':
        data[:, 96:] = np.nan

    if name == 'ILTJ130634.72+553657.7':
        data[107:, :] = np.nan

    if name == 'ILTJ130804.04+550835.4':
        data[:, :46] = np.nan
        data[:, 70:] = np.nan
        data[:13, :] = np.nan

    if name == 'ILTJ130821.80+562744.1':
        data[:, :23] = np.nan
        data[:, 100:] = np.nan

    if name == 'ILTJ130857.14+542915.9':
        data[:, 98:] = np.nan
        data[:, :23] = np.nan

    if name == 'ILTJ130907.51+545520.4':
        data[:, 96:] = np.nan

    if name == 'ILTJ130926.29+534820.0':
        data[54:57, 76:79] = np.nan
        data[44:53, 69:74] = np.nan

    if name == 'ILTJ131115.53+534356.8':
        data[101:, :] = np.nan

    if name == 'ILTJ131236.04+535524.7':
        data[:, 83:] = np.nan

    if name == 'ILTJ125757.12+554224.1':
        data[:, :40] = np.nan
        data[:, 71:] = np.nan

    if name == 'ILTJ125840.51+534935.5':
        data[79:, :] = np.nan
        # d[:45, :] = np.nan
        data[:, 70:] = np.nan

    if name == 'ILTJ125852.58+555527.3':
        data[:, :21] = np.nan
        data[71:, :] = np.nan

    if name == 'ILTJ130005.50+551315.1':
        data[:, :32] = np.nan

    if name == 'ILTJ130112.93+550410.5':
        data[:33, :] = np.nan

    if name == 'ILTJ130148.82+544723.9':
        data[:50, :] = np.nan

    if name == 'ILTJ130331.44+534400.9':
        data[:, :43] = np.nan
        data[67:85, 63:88] = np.nan
        data[71:, :] = np.nan
        data[:, 72:] = np.nan
        data[65:85, 66:88] = np.nan

    if name == 'ILTJ130441.11+551952.7':
        data[:, :39] = np.nan
        data[:, 73:] = np.nan
        data[:37, :] = np.nan

    if name == 'ILTJ130531.20+535432.3':
        data[:, :39] = np.nan
        data[:, 97:] = np.nan
        data[:30, :] = np.nan
        data[:30, :] = np.nan
        data[45:53, 55:64] = np.nan

    if name == 'ILTJ130605.63+555127.6':
        data[:, :41] = np.nan

    if name == 'ILTJ130634.72+553657.7':
        data[73:85, 25:42] = np.nan
        data[30:41, 43:56] = np.nan

    if name == 'ILTJ130638.74+541448.8':
        data[:, :22] = np.nan
        data[:, 91:] = np.nan
        data[:45, :] = np.nan

    if name == 'ILTJ130749.04+545120.0':
        data[:, :11] = np.nan
        data[91:, :] = np.nan
        data[:42, :] = np.nan
        data[42:51, 54:63] = np.nan

    if name == 'ILTJ130804.04+550836.4':
        data[:17, :] = np.nan

    if name == 'ILTJ130857.14+542915.9':
        data[58:68, 35:45] = np.nan
        data[66:77, 35:52] = np.nan
        data[:, 71:] = np.nan
        data[40:49, 60:76] = np.nan

    if name == 'ILTJ130926.29+534820.0':
        data[65:77, 42:55] = np.nan

    if name == 'ILTJ130331.30+540024.1':
        data[174:, :] = np.nan
        data[:, :50] = np.nan

    if name == 'ILTJ104416.22+524037.3':
        data[:, 120:] = np.nan

    if name == 'ILTJ104438.12+521539.4':
        data[100:, :] = np.nan

    if name == 'ILTJ104611.79+522931.2':
        data[:, :130] = np.nan
        data[:, 170:] = np.nan

    if name == 'ILTJ104647.32+531352.7':
        data[:, 100:] = np.nan

    if name == 'ILTJ104650.69+544807.2':
        data[:, :70] = np.nan
        data[:, 130:] = np.nan

    if name == 'ILTJ104703.03+523023.4':
        data[:20, :] = np.nan
        data[70:, :] = np.nan

    if name == 'ILTJ104713.87+530240.3':
        data[170:, :] = np.nan
        data[:100, :] = np.nan

    if name == 'ILTJ104725.20+480043.4':
        data[:, :20] = np.nan

    if name == 'ILTJ104733.68+530653.6':
        data[:60, :] = np.nan
        data[120:, :] = np.nan

    if name == 'ILTJ104747.59+521414.9':
        data[:, 77:] = np.nan
        data[100:, :] = np.nan

    if name == 'ILTJ104843.46+482830.4':
        data[:, 48:] = np.nan

    if name == 'ILTJ104948.06+520821.8':
        data[:40, :] = np.nan

    if name == 'ILTJ104611.79+522931.2':
        data[:, :130] = np.nan
        data[:, 170:] = np.nan

    if name == 'ILTJ104706.09+534417.1':
        data[:88, 100:] = np.nan

    if name == 'ILTJ104713.87+530240.3':
        data[:, 250:] = np.nan

    if name == 'ILTJ104747.59+521414.9':
        data[100:, :] = np.nan

    if name == 'ILTJ104955.05+461923.0':
        data[:, :20] = np.nan

    if name == 'ILTJ105028.76+544849.7':
        data[30:, :] = np.nan

    if name == 'ILTJ105043.48+471027.2':
        data[:, :60] = np.nan
        data[100:, :] = np.nan
        data[:60, :] = np.nan

    if name == 'ILTJ105045.37+463949.2':
        data[80:, :] = np.nan

    if name == 'ILTJ105047.02+523523.8':
        data[:, :45] = np.nan

    if name == 'ILTJ105106.79+482535.0':
        data[:24, :] = np.nan

    if name == 'ILTJ105122.62+473619.6':
        data[:, 100:] = np.nan

    if name == 'ILTJ105144.36+515814.4':
        data[:, 100:] = np.nan

    if name == 'ILTJ105224.33+524345.5':
        data[:, :16] = np.nan
        data[60:, :] = np.nan

    if name == 'ILTJ105238.10+484019.7':
        data[80:, :] = np.nan

    if name == 'ILTJ105239.07+463413.9':
        data[66:, 75:] = np.nan

    if name == 'ILTJ105252.33+553845.6':
        data[:, :30] = np.nan

    if name == 'ILTJ105329.78+484129.3':
        data[:30, :] = np.nan
        data[80:, :] = np.nan
        data[:, 100:] = np.nan

    if name == 'ILTJ105351.79531021.2':
        data[70:, :] = np.nan
        data[:34, 54:] = np.nan

    if name == 'ILTJ105418.83+512447.7':
        data[:20, :] = np.nan
        data[:, :32] = np.nan

    if name == 'ILTJ105430.25+542506.2':
        data[:, 140:] = np.nan

    if name == 'ILTJ105444.97+561718.5':
        data[:45, :] = np.nan

    if name == 'ILTJ105456.44+533546.6':
        data[:, 49:] = np.nan

    if name == 'ILTJ105510.77+461609.0':
        data[:, 120:] = np.nan
        data[137:, :] = np.nan
        data[:40, :] = np.nan

    if name == 'ILTJ105514.67+490001.8':
        data[60:, :] = np.nan
        data[:, :26] = np.nan

    return data


def measure_extent_radio_galaxies(
        my_dir='/data5/sean/deep-fields/catalogues/', nah=False):
    """Measure the size of radio galaxies.
    """

    results_csv = f'{my_dir}../results/extention-radio-galaxies.csv'
    if nah:
        return results_csv  # do not do anything

    df = pd.read_csv(f'{my_dir}radio.galaxies.ldr1.csv')  # load data
    fields = [f'/data1/lotss-data/{m}-mosaic.fits' for m in df['Mosaic_ID']]
    thresholds = 5 * df['Isl_rms'] / 1000  # converting to jansky
    plt.figure(figsize=(8, 8))

    result_header = ('Name,Class,LM flux (mJy),LM size ("),Width ("),Width' +
                     ' (kpc),5 * rms (mJy),Redshift,NAT,WAT,D-D\n')

    if not os.path.exists(results_csv):
        with open(results_csv, 'a') as f:
            f.write(result_header)

    # skipping these sources as they are too complicated
    skip_list = ['ILTJ105155.09+552329.0', 'ILTJ104703.03+523023.4']
    i = 0

    for (source_name, ra, dec, field, threshold, fri, frii, redshift, nat, wat,
         dd, lm_size, lm_flux) in zip(df['Source_Name_1'], df['RA_1'],
                                      df['DEC_1'], fields, thresholds,
                                      df['FR1'], df['FR2'], df['z_best'],
                                      df['NAT'], df['WAT'], df['D-D'],
                                      df['LM_dec_size'], df['LM_Flux']):

        save = f'{my_dir}../images/extention-radio-galaxies/{source_name}.png'
        if os.path.exists(save):
            print(f'{save} already exists so it is being skipped.')
            continue

        if source_name in skip_list:
            print(f'{source_name} is too complicated so it is being skipped.')
            continue

        source_type = 'FR-I' if fri else 'FR-II'  # already removed small
        # sources and ambiguous sources from the sample

        # if source_name != 'ILTJ130804.04+550835.4':  # for testing one source
        #     continue

        # try:
        hdu = fits.open(field)[0]
        # except FileNotFoundError:  # fits might have different name
        #     print(f'{field} does not exist.')
        #     continue

        wcs = WCS(hdu.header, naxis=2)
        sky_position = SkyCoord(ra, dec, unit='deg')
        s = math.ceil((lm_size * 2) / 60)
        # already have adhoc masks for P196+55 field at [3, 3] arcmin but for
        # the rest we take the size of the source from lomorph
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

        # this removes highly curved sources
        # if lm_size / asec_max > 1.75 or asec_max / lm_size > 1.75:
        #     print(f'{source_name} is too curved so it is being skipped.')
        #     continue

        # for good cells find the distance from it to the point on the line, we
        # assume we look down the longest line and the source is flat in the
        # plane of the sky (https://stackoverflow.com/a/52756183/6386612)
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
        asec_width = (my_max_width + np.abs(my_min_width)) * 1.5  # 1.5"/pixel

        plt.imshow(d, vmin=0, vmax=np.nanmax(d), origin='lower',
                   norm=DS9Normalize(stretch='arcsinh'))

        plt.colorbar()
        plt.plot([max_y1, max_y2], [max_x1, max_x2], color='red', alpha=0.5,
                 lw=2)
        plt.plot([width_y_max], [width_x_max], marker='o', markersize=10,
                 color='red', alpha=0.5)
        plt.plot([width_y_min], [width_x_min], marker='o', markersize=10,
                 color='red', alpha=0.5)

        plt.minorticks_on()
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

        _, kpc = get_dl_and_kpc_per_asec(z=redshift)
        kpc_width = kpc * asec_width

        results = (f'{source_name},{source_type},{lm_flux},{lm_size},' +
                   f'{asec_width},{kpc_width},{threshold * 1000},{redshift},' +
                   f'{nat},{wat},{dd}\n')

        with open(results_csv, 'a') as f:
            f.write(results)

        i += 1
        if i > 100:
            import sys
            sys.exit()
        # build clean sample on visual inspection


def main():
    measure_extent_radio_galaxies(nah=False)


if __name__ == '__main__':
    main()
