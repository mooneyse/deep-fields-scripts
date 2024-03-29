#!/usr/bin/env python3

"""Fit a Gaussian point spread function to a point source and subtract it from
a source with diffuse emission.
"""

import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.table import Table
from astropy.wcs import WCS
from scipy.ndimage.interpolation import map_coordinates, shift
import scipy.optimize as opt
from ds9norm import DS9Normalize

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '10 October 2019'


def get_df(filename, format, index):
    """Get the dataframe.

    Parameters
    ----------
    filename : string
        Name of the file.
    format : string
        The format of the file.
    index : integer
        The index.

    Returns
    -------
    Pandas dataframe
        Contents of the file.
    """
    if format == 'csv':
        df = pd.read_csv(filename)
    else:
        data = Table.read(filename, format=format)
        df = data.to_pandas()
    df.set_index(index, inplace=True)
    return df


def get_position(df, cat_dir):
    """Look up the position of the blazar.

    Parameters
    ----------
    df : Pandas dataframe
        Name of the file.
    cat_dir : string
        Directory of the catalogue.

    Returns
    -------
    list
        Source names.
    list
        Positions (RA, declination).
    list
        Catalogues.
    list
        FITS images.
    """
    # df = df[df['Field'] == field]  # filters to one field
    blazar_names = df.index.tolist()
    blazar_positions, catalogues, fits_images = [], [], []
    for ra, dec, field in zip(df['RA'], df['DEC'], df['Field']):
        blazar_positions.append([ra, dec])
        field = field.lower().replace(' ', '.')
        catalogues.append(f'{cat_dir}/{field}.11.06.2019.cat.fits')
        fits_images.append(f'{cat_dir}/{field}.11.06.2019.img.fits')
    return blazar_names, blazar_positions, catalogues, fits_images


def nearest_point_source(df, position, s_code='S', flux_threshold=0.01,
                         distance_threshold=1.5, elongation_threshold=1.2):
    """Find the nearest bright point source to a given blazar.

    Parameters
    ----------
    df : Pandas dataframe
        The data.
    position : list
        Position of the source.
    s_code : string, optional
        PyBDSF source code to use. The default is 'S'.
    flux_threshold : float, optional
        Minimum flux of a source to accept, in Jy per beam to use, when looking
        for a nearby source. The default is 10 mJy.
    distance_threshold : float, optional
        Maximum distance the source can be to the target, in degrees. The
        default is 1.5.
    elongation_threshold : float, optional
        How elongated the point source can be and still accepted. If a source
        is too long, it may be intrinsic, as the beam is relatively
        symmetrical. It is given as a ratio betwen the major and minor axes.

    Returns
    -------
    float
        Distance to the target source.
    list
        Position (RA, declination) of the point source.
    """
    # disable warning, see https://stackoverflow.com/a/20627316/6386612
    pd.options.mode.chained_assignment = None
    df.reset_index(inplace=True)  # remove Source_id as index column
    df['S_Code'] = df['S_Code'].str.decode('utf-8')  # byte string encoded
    df_point_sources = df[(df['S_Code'] == s_code) &
                          (df['Total_flux'] > flux_threshold)]
    blazar = SkyCoord(position[0], position[1], unit='deg')
    separations, elongations = [], []
    for ra, dec, major, minor in zip(df_point_sources['RA'],
                                     df_point_sources['DEC'],
                                     df_point_sources['Maj'],
                                     df_point_sources['Min']):
        point_source = SkyCoord(ra, dec, unit='deg')
        separation = blazar.separation(point_source)
        separations.append(separation.deg)
        elongations.append(major / minor)
    df_point_sources['Separation'] = separations
    df_point_sources['Elongation'] = elongations
    df_point_sources = df_point_sources[df_point_sources['Elongation'] <
                                        elongation_threshold]
    nearest = df_point_sources.loc[df_point_sources['Separation'].idxmin()]
    point_source_position = [nearest['RA'], nearest['DEC']]
    results = {'i': 'ILTJ' + str(nearest['Source_id']),
               's': nearest['Separation'], 'f': nearest['Total_flux'] * 1000}
    print('{i} is {s:.2f} degrees away and has a total flux density of ' +
          '{f:.2f} mJy.'.format(**results))
    if results['s'] > distance_threshold:
        print('The point source is too far away.')
        return False, False
    return results['i'], point_source_position


def get_fits(filename):
    """Open FITS image.

    Parameters
    ----------
    filename : string
        Name of the file.

    Returns
    -------
    Astroy object
        HDU.
    Astropy object.
        World coordinate system.
    """
    hdu = fits.open(filename)[0]
    wcs = WCS(hdu.header, naxis=2)
    return hdu, wcs


def get_data(position, hdu, wcs, size=[2, 2] * u.arcmin):
    """Cut out the source from the FITS data.

    Parameters
    ----------
    position : list
        Position of the source, as RA and declination in degrees.
    hdu : Astropy object
        The FITS HDU.
    wcs : Astropy object
        The FITS world coordinate system.
    size : list, optional
        The size of the desired cut out. The default is 2 by 2 arcminutes.

    Returns
    -------
    pandas dataframe
        Contents of the file.
    """
    sky_position = SkyCoord(position[0], position[1], unit='deg')
    cutout = Cutout2D(np.squeeze(hdu.data), sky_position, size=size, wcs=wcs)
    data = cutout.data
    return data


def housekeeping(name, data):
    '''Some ad hoc adjustments to a few sources, made after visual inspection.'''
    if name == '5BZQJ1422+3223':
        print('Doing a little housekeeping on {}.'.format(name))
        data[42:47, 55:60] = 0
        data[6:10, 74:79] = 0

    if name == '5BZBJ1426+3404':
        print('Doing a little housekeeping on {}.'.format(name))
        data[20:28, 23:32] = 0  # data[:10, :10] = 0
        data[68:75, 33:42] = 0

    if name == '5BZQJ1429+3529':
        print('Doing a little housekeeping on {}.'.format(name))
        data[40:48, 10:17] = 0
        data[10:26, 40:60] = 0  # data[:10, 20:] = 0

    if name == '5BZQJ1435+3353':
        print('Doing a little housekeeping on {}.'.format(name))
        data[38:45, 5:11] = 0  # data[8:13, 28:33] = 0
        data[73:79, 54:60] = 0
        data[27:33, 48:53] = 0
        data[16:22, 64:71] = 0

    if name == '5BZQJ1437+3618':
        print('Doing a little housekeeping on {}.'.format(name))
        data[0:31, 4:26] = 0  # data[1:10, :3] = 0

    if name == '5BZQJ1437+3519':
        print('Doing a little housekeeping on {}.'.format(name))
        # data[:, 60:] = 0

    if name == '5BZBJ1558+5625':
        print('Doing a little housekeeping on {}.'.format(name))
        data[65:70, 36:43] = 0

    if name == '5BZBJ1605+5421':
        print('Doing a little housekeeping on {}.'.format(name))
        data[53:64, 16:27] = 0
        data[72:79, 35:45] = 0
        data[76:79, 76:80] = 0
        data[48:56, 75:80] = 0
        data[6:11, 70:76] = 0
        data[24:31, 42:48] = 0

    if name == '5BZQJ1606+5405':
        print('Doing a little housekeeping on {}.'.format(name))
        data[57:68, 69:79] = 0
        data[25:30, 11:17] = 0

    if name == '5BZQJ1608+5613':
        print('Doing a little housekeeping on {}.'.format(name))
        data[19:30, 21:37] = 0

    if name == '5BZQJ1619+5256':
        print('Doing a little housekeeping on {}.'.format(name))
        data[42:46, 48:55] = 0
        data[27:32, 35:42] = 0

    if name == '5BZBJ1037+5711':
        print('Doing a little housekeeping on {}.'.format(name))
        data[21:26, 27:32] = 0
        data[55:61, 39:46] = 0

    return data


def regrid(data, new_size=10, normalise=True):
    '''Map the data onto a larger array.'''
    l, w = data.shape
    new_l = l * new_size
    new_w = w * new_size
    new_dimensions = []
    for old_shape, new_shape in zip(data.shape, (new_l, new_w)):
        new_dimensions.append(np.linspace(0, old_shape - 1, new_shape))
    coordinates = np.meshgrid(*new_dimensions, indexing='ij')
    new_data = map_coordinates(data, coordinates)
    new_data = new_data / np.max(new_data) if normalise else new_data
    return new_data  # this function was tested and worked as expected


def gaussian(xy, amplitude, x0, y0, sigma_x, sigma_y, theta, offset):
    '''Fit a two-dimensional Gaussian to the data.'''
    x, y = xy
    a = ((np.cos(theta) ** 2) / (2 * sigma_x ** 2) + (np.sin(theta) ** 2) / (2 * sigma_y ** 2))
    b = (-(np.sin(2 * theta)) / (4 * sigma_x ** 2) + (np.sin(2 * theta)) / (4 * sigma_y ** 2))
    c = ((np.sin(theta) ** 2)/(2 * sigma_x ** 2) + (np.cos(theta) ** 2) / (2 * sigma_y ** 2))
    g = (offset + amplitude * np.exp( -(a * ((x - x0) ** 2) + 2 * b * (x - x0) * (y - y0) + c * ((y - y0) **2 ))))
    g_ravel = g.ravel()
    return g_ravel


def make_model(data, amplitude=1, sigma_x=30, sigma_y=30, theta=0, offset=0):
    '''Fit a model to the data.'''
    len_x, len_y = data.shape
    range_x = range(len_x)
    range_y = range(len_y)
    x, y = np.meshgrid(range_x, range_y)
    x0 = len_x / 2
    y0 = len_y / 2
    data_ravel = data.ravel()
    p0 = (amplitude, x0, y0, sigma_x, sigma_y, theta, offset)
    popt, pcov = opt.curve_fit(gaussian, (x, y), data_ravel, p0=p0)
    model = gaussian((x, y), *popt)
    model_reshape = model.reshape(len_x, len_y)
    return model_reshape


def match_peaks(data, model, cval=0):
    '''Shift the data so that the index of the maximum value in the data matches up with that of the model. This is needed to ensure accurate subtraction.'''
    data_max = data.argmax()
    model_max = model.argmax()
    data_shape = data.shape
    model_shape = model.shape
    data_peak_row, data_peak_column = np.unravel_index(data_max, data_shape)
    model_peak_row, model_peak_column = np.unravel_index(model_max, model_shape)
    shifting = (model_peak_row - data_peak_row, model_peak_column - data_peak_column)
    shifted_data = shift(data, shift=shifting, cval=cval)
    numbers = {'d': (data_peak_row, data_peak_column), 's': shifting, 'm': (model_peak_row, model_peak_column)}
    print('The data were shifted from {d} by {s} to {m}.'.format(**numbers))
    return shifted_data


def get_noise_catalogue(df, blazar_name, new_data='', sigma=5):
    '''Calculate the noise in the blazar image.'''
    rms = df.loc[blazar_name, 'Isl_rms']
    five_sigma = [rms * sigma]
    return five_sigma


def jy_per_beam_to_jy(flux, beam=6, sigma=5):
    '''Convert Jansky per beam measurements to Jansky measurements.'''
    bmaj = beam * u.arcsec
    bmin = beam * u.arcsec
    fwhm_to_sigma = 1 / np.sqrt(8 * np.log(2))
    beam_area = 2 * np.pi * bmaj * bmin * (fwhm_to_sigma ** 2)
    beam_area = beam_area.to(u.sr)
    flux_jy_per_beam = flux * u.Jy / u.sr
    flux_jy = flux_jy_per_beam * beam_area
    # BUG this is way off but the code looks right
    # see http://docs.astropy.org/en/stable/api/astropy.units.equivalencies.brightness_temperature.html
    # also see https://astronomy.stackexchange.com/a/20391/20014
    # also see my note from 2019-06-20
    return flux_jy #


def diffuse_fraction(df, name, blazar, diffuse, threshold):
    '''Calculate what the fraction of emission is diffuse.'''
    blazar_catalogue = df.loc[name, 'Total_flux']
    blazar_image = blazar[blazar > threshold].sum()
    # print(f'There is {blazar_image} Jy per beam above the given \N{GREEK SMALL LETTER SIGMA}.')
    diffuse_image = diffuse[diffuse > threshold].sum()
    fraction = diffuse_image / blazar_image
    blazar_diffuse_catalogue = blazar_catalogue * fraction
    blazar_core_catalogue = blazar_catalogue * (1 - fraction)
    print('{0:.2%} of the emission is not from the core.'.format(fraction))
    results = {'total': blazar_catalogue * 1000, 'core': blazar_core_catalogue * 1000, 'diffuse': blazar_diffuse_catalogue * 1000}  # mJy
    print('The total flux density is {total:.2f} mJy, with the core contributing {core:.2f} mJy and the diffuse emission contributing {diffuse:.2f} mJy.'.format(**results))
    return name, blazar_diffuse_catalogue


def make_plot(position, data, title, rows=2, columns=7, origin='lower',
              vmin=0, vmax=1, zmin=0, zmax=1, axis='off', elev=24, azim=29,
              linewidth=0.25, pane_colour='white', line_colour='black',
              cmap='magma_r', projection='3d', stretch='arcsinh', levels='',
              plot='', layer='', contour_colours=['blue', 'magenta'], pad=0.05,
              size='5%', orientation='horizontal', location='bottom',
              linestyles='dashed'):
    '''Plot the image on a grid.'''
    ax0 = plt.subplot(rows, columns, position)
    ax0.set_title(title)
    ax0.get_xaxis().set_ticks([])
    ax0.get_yaxis().set_ticks([])
    image_size = 60  # image is 1 * 1 arcminutes
    len_x, len_y = data.shape
    x0, x1 = len_x * 0.65, len_x * 0.9
    y = len_y * 0.1
    plt.text((x0 + x1) / 2, y * (1.4), str(int(image_size * 0.25)) + '"', horizontalalignment='center', verticalalignment='center', fontsize=12)
    plt.plot([x0, x1], [y, y], 'k-', lw=1.5)
    image = ax0.imshow(data, origin=origin, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=DS9Normalize(stretch=stretch))
    divider = make_axes_locatable(ax0)
    cax0 = divider.append_axes(location, size=size, pad=pad)
    plt.colorbar(image, cax=cax0, orientation=orientation)
    if plot is 'blazar':
        ax0.contour(data, levels=levels, origin=origin, colors=contour_colours[0])
    elif plot is 'diffuse':
        ax0.contour(layer, levels=levels, origin=origin, colors=contour_colours[0])
        ax0.contour(data, levels=levels, origin=origin, linestyles=linestyles, colors=contour_colours[1])
    len_x, len_y = data.shape
    range_x = range(len_x)
    range_y = range(len_y)
    ax1 = plt.subplot(rows, columns, int(position + columns), projection=projection)
    x, y = np.meshgrid(range_x, range_y)
    ax1.get_xaxis().set_ticks([])
    ax1.get_yaxis().set_ticks([])
    ax1.zaxis.set_tick_params(pad=10)  # ax1.set_zticks([]) to remove z labels
    ax1.xaxis.pane.fill = False
    ax1.yaxis.pane.fill = False
    ax1.zaxis.pane.fill = False
    ax1.xaxis.pane.set_edgecolor(pane_colour)
    ax1.yaxis.pane.set_edgecolor(pane_colour)
    ax1.zaxis.pane.set_edgecolor(pane_colour)
    ax1.grid(False)
    ax1.set_zlim(zmin, vmax)
    ax1.view_init(elev=elev, azim=azim)
    ax1.plot_surface(x, y, data, cmap=cmap, vmin=vmin, vmax=vmax, linewidth=linewidth, edgecolors=line_colour)


def new_plots(pos, data, layer=None, plot_type='blazar', vmax=1, levels=[]):
    ax0 = plt.subplot(1, 2, pos)
    ax0.get_xaxis().set_ticks([])
    ax0.get_yaxis().set_ticks([])
    image_size = 60 * 2  # image is 2 * 2 arcminutes
    len_x, len_y = data.shape
    x0, x1 = len_x * (0.65 + 0.125), len_x * 0.9
    y = len_y * 0.1
    plt.text((x0 + x1) / 2, y * (0.7), str(int(image_size * 0.25 * 0.5)) + '"', horizontalalignment='center', verticalalignment='center', fontsize=15, color='w')
    plt.plot([x0, x1], [y, y], 'w-', lw=1.5)

    image = ax0.imshow(data, origin='lower', cmap='viridis', vmax=vmax, vmin=0, norm=DS9Normalize(stretch='arcsinh'))
    circle1 = plt.Circle([8, 8], (x1 - x0) * (6 / 15) / 2, ls='--', fc='None', ec='w')
    ax0.add_artist(circle1)
    divider = make_axes_locatable(ax0)
    cax0 = divider.append_axes('bottom', size='5%', pad=0.05)
    cbar = plt.colorbar(image, cax=cax0, orientation='horizontal')
    cbar.ax.tick_params(labelsize=15)
    cbar.ax.set_xlabel(r'Jy beam$^{-1}$', fontsize=15)
    if plot_type is 'blazar':
        ax0.contour(data, levels=levels, origin='lower', colors='white', alpha=0.67)
    elif plot_type is 'diffuse':
        ax0.contour(layer, levels=levels, origin='lower', colors='white', alpha=0.67)
        ax0.contour(data, levels=levels, origin='lower', colors='magenta')  # in this case, data is the core subtracted array
        # ax0.contourf(data, levels=[levels[0],levels[0]*100], origin='lower', colors='magenta', alpha=0.2)


def main():
    '''Fit a Gaussian point spread function to a point source and subtract it from a source with diffuse emission.'''
    formatter_class = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=formatter_class)
    # parser.add_argument('-i', '--image', required=False, type=str, default='/mnt/closet/ldr2-blazars/deep-fields/bootes-image.fits', help='FITS image of the field')
    parser.add_argument('-c', '--catalogue_dir', required=False, type=str, default='/mnt/closet/deep-fields/catalogues', help='Directory with the catalogues')
    parser.add_argument('-d', '--csv', required=False, type=str, default='/mnt/closet/deep-fields/catalogues/deep.fields.11.06.2019.cat.csv', help='CSV catalogue of the blazars')
    parser.add_argument('-o', '--output', required=False, type=str, default='/mnt/closet/deep-fields/images/core-subtract', help='Directory to save the plots')

    args = parser.parse_args()
    # image = args.image
    catalogue_dir = args.catalogue_dir
    csv = args.csv
    output = args.output

    font = 'STIXGeneral'
    math_font = 'cm'
    font_size = 12
    figsize = (17, 8)  # (40, 20)
    bbox_inches = 'tight'
    testing = False
    new_size = 10
    my_blazars, my_diffuse = [], []

    df_blazars = get_df(csv, format='csv', index='Source name')
    blazar_names, blazar_positions, catalogues, images = get_position(df_blazars, cat_dir=catalogue_dir)
    for i, (blazar_name, blazar_position, catalogue, image) in enumerate(zip(blazar_names, blazar_positions, catalogues, images)):
        if testing:
            if i != 0:  # do one at a time
                sys.exit()
        # if blazar_name != '5BZQJ1437+3519':
        #     continue
        print(f'Analysing {blazar_name} which is in {image} (blazar {i + 1} of {len(blazar_names)}).')
        df_cat = get_df(catalogue, format='fits', index='Source_id')
        point_source_id, point_source_position = nearest_point_source(df_cat, blazar_position)
        if point_source_id is False:
            print('Going to the next iteration.')
            continue
        hdu, wcs = get_fits(filename=image)
        blazar_data = get_data(position=blazar_position, hdu=hdu, wcs=wcs)
        blazar_data = housekeeping(blazar_name, blazar_data)
        point_source_data = get_data(position=point_source_position, hdu=hdu, wcs=wcs)
        blazar_regrid = regrid(blazar_data, new_size=new_size, normalise=False)  # peak and total values change with regridding
        point_source_regrid = regrid(point_source_data, new_size=new_size, normalise=False)
        model = make_model(point_source_regrid, sigma_x=4, sigma_y=4)
        point_source_residual = point_source_regrid - model
        scaled_model = (model * np.max(blazar_regrid) / np.max(point_source_regrid))
        blazar_shifted = match_peaks(blazar_regrid, scaled_model)
        blazar_residual = blazar_shifted - scaled_model
        blazar_regrid_back = regrid(blazar_shifted, new_size=1 / new_size, normalise=False)  # regrid the blazar and blazar residual data back to the native resolution
        blazar_residual_regrid_back = regrid(blazar_residual, new_size=1 / new_size, normalise=False)
        five_sigma = get_noise_catalogue(df_blazars, blazar_name)
        if blazar_name == '5BZQJ1437+3519':
            print('Doing a little more housekeeping on {}.'.format(blazar_name))
            blazar_regrid_back[:, 60:] = 0
            blazar_residual_regrid_back[:, 60:] = 0
        b, d = diffuse_fraction(df=df_blazars, name=blazar_name, blazar=blazar_regrid_back, diffuse=blazar_residual_regrid_back, threshold=five_sigma)
        my_blazars.append(b)
        my_diffuse.append(d)
        continue  # skip the plotting as I have that already
        savefig = output + '/' + blazar_name + '.png'
        matplotlib.rcParams['font.family'] = font
        matplotlib.rcParams['mathtext.fontset'] = math_font
        matplotlib.rcParams['font.size'] = font_size
        plt.figure(figsize=figsize)
        # make_plot(position=1, data=point_source_regrid, title=point_source_id, vmax=np.max(point_source_regrid))
        # make_plot(position=2, data=model, title=point_source_id + ' model', vmax=np.max(point_source_regrid))
        # make_plot(position=3, data=point_source_residual, title=point_source_id + ' residual', vmax=np.max(point_source_regrid))
        # make_plot(position=4, data=blazar_shifted, title=blazar_name, vmax=np.max(blazar_shifted))
        # make_plot(position=5, data=blazar_residual, title=blazar_name + ' diffuse', vmax=np.max(blazar_shifted))
        # make_plot(position=6, data=blazar_regrid_back, title=blazar_name, levels=five_sigma, plot='blazar', vmax=np.max(blazar_regrid_back))
        # make_plot(position=7, data=blazar_residual_regrid_back, title=blazar_name + ' diffuse', levels=five_sigma, plot='diffuse', layer=blazar_regrid_back, vmax=np.max(blazar_regrid_back))
        # if testing:
        #     plt.show()
        # else:
        #     plt.savefig(savefig, bbox_inches=bbox_inches)
        #     print(f'Done! The plot is saved. View it with this: gpicview {savefig}')
        # TODO make a df from the results and export it as a csv
        # TODO make a new plotter only showing the 2d blazar_regrid_back and blazar_residual_regrid_back with viridis
        #      could I get this contour onto the radio image with proper axes?
        #      my results might not be correct as I am manually blocking out areas for certain sources
        # TODO plot negative flux to see if there are any bowls, using a diverging colour scale
        # TODO does the flux within 5 sigma equal the catalogue flux?
        new_plots(pos=1, data=blazar_regrid_back, vmax=np.max(blazar_regrid_back), levels=five_sigma)
        # old_max = np.max(blazar_regrid_back)
        # blazar_regrid_back[24:43, 36:42] = 0
        # blazar_regrid_back[38:40, 42:43] = 0
        # blazar_regrid_back[35:38, 35:36] = 0
        # blazar_regrid_back[27:30, 42:43] = 0
        # new_plots(pos=2, data=blazar_regrid_back, vmax=old_max, levels=five_sigma, layer=blazar_regrid_back, plot_type='diffuse')
        new_plots(pos=2, data=blazar_residual_regrid_back, vmax=np.max(blazar_regrid_back), levels=five_sigma, layer=blazar_regrid_back, plot_type='diffuse')
        # plt.title(blazar_name, color='white')
        plt.tight_layout()
        plt.show()
        # plt.savefig(f'{output}/2-plots-{blazar_name}.png')

    print()
    for b, d in zip(my_blazars, my_diffuse):
        print(f'{b} {d}')


if __name__ == '__main__':
    main()
