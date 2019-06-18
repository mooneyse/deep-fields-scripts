#!/usr/bin/env python3

'''Plot the Deep Fields radio contours over the SDSS optical images which are
fetched from the website with an API.'''

import warnings
warnings.filterwarnings('ignore')  # supress warnings

import argparse
import aplpy
import glob
import sys
import pandas as pd
import matplotlib.pyplot as plt
from astroquery.sdss import SDSS
from astropy import coordinates as coords
from pathlib import Path

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '17 June 2019'

def query_sdss(ra, dec, unit='deg', frame='icrs', band='i', spectro=False):
    '''Query SDSS. Many filters are available. See this URL for more:
    http://skyserver.sdss.org/dr2/en/proj/advanced/color/sdssfilters.asp.'''

    position = coords.SkyCoord(ra, dec, unit=unit, frame=frame)
    id = SDSS.query_region(position, spectro=spectro)
    images = SDSS.get_images(matches=id, band=band)

    return images[0]  # use first image, which is the corrected one I think


def make_plot(fits_file, df, format='png', north=True, figsize=12, radius=120,
              cmap='bone_r', vmin=0, radio_cmap='YlOrRd_r', color='white',
              radio_levels=[4, 8, 16, 32], length=1074, width=1009, vmax=1,
              save_directory='/mnt/closet/ldr2-blazars/images/contour',
              fill='black'):
    '''Plot contours from one FITS file over the image of another. Radius is
    given in arcseconds.'''

    # house-keeping
    blazar_name = fits_file.split('/')[-1].split('-P')[0]
    blazar_row = df[df['name'] == blazar_name]
    ra, dec = float(blazar_row['ra']), float(blazar_row['dec'])
    save = '{}/{}.{}'.format(save_directory, blazar_name, format)
    sdss_data = True
    if len(glob.glob(save_directory + '/*')) > 300:
        print('300 files reached!')
        return

    if Path(save).is_file():  # allows the script to be restarted
        print('{} already exists so it is being skipped.'.format(save))
        return

    try:
        image_file = query_sdss(ra, dec)
        image = aplpy.FITSFigure(image_file, north=north, figsize=(figsize, figsize))
    except:  # if there is no sdss just plot the radio data
        print('SDSS match for {} not found.'.format(blazar_name))
        sdss_data = False
        image = aplpy.FITSFigure(fits_file, figsize=(figsize, figsize))

    # make image
    image.show_colorscale(cmap=cmap, vmin=vmin, vmax=vmax)
    image.recenter(ra, dec, radius=radius / (60 * 60))  # degrees
    image.add_scalebar(radius / (60 * 60 * 2))
    image.scalebar.set_label(str(int(radius / 2)) + '"')
    # image.show_grid()
    image.set_title(blazar_name if sdss_data else blazar_name + ' (SDSS data not found)')
    image.show_contour(fits_file, cmap=radio_cmap, levels=radio_levels)  # mJy
    image.save(save)
    print('Image saved at {}.'.format(save))


def radio_over_optical(radio_directory, bzcat_csv):
    '''Get the list of sources and for each, call a function to plot the radio
    contours over the optical image.'''

    df = pd.read_csv(bzcat_csv)  # get catalogue which has positions we need
    df.rename(columns={' Source name ': 'name', ' RA (J2000.0) ': 'ra',
                       ' Dec (J2000.0) ': 'dec'}, inplace=True)
    df['name'] = df['name'].str.strip()  # remove spaces from blazar names
    fits_files = glob.glob('{}/*'.format(radio_directory))  # get fits files

    for fits_file in fits_files:
        make_plot(fits_file, df)


def do_plotting(name, ra, dec, noise, optical, radio):
    '''Given FITS files with radio and optical data, this will plot the radio
    contours over the optical image.

    Args:
        name (str): The name of the source, used for saving the file.
        ra (float): The source right ascension.
        dec (float): The source declination.
        noise (float): The RMS noise to be used to calculate the contours.
        optical (str): FITS file containing data for the map.
        radio (str): FITS file containing data for the contours.

    Returns:
        NoneType: Returns nothing.'''

    cmap = 'Greys'
    radius = 1/60
    contours = [5, 10, 20, 40, 80, 160, 320, 640, 1280]
    colors = ['#e6194b', '#f58231', '#ffe119', '#bfef45', '#3cb44b', '#42d4f4', '#4363d8', '#911eb4', '#f032e6']

    image = aplpy.FITSFigure(optical, north=True)  # figsize=(12, 12)
    image.recenter(ra, dec, radius=radius)

    image.show_colorscale(cmap=cmap, vmin=0, vmax=1)
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['mathtext.fontset'] = 'dejavuserif'

    image.add_beam(6 * 1/60/60, 6 * 1/60/60, 0)
    image.beam.set_linestyle('dashed')
    image.beam.set_edgecolor('black')
    image.beam.set_facecolor('none')

    levels = [contour * noise for contour in contours]
    image.show_contour(radio, levels=levels, colors=colors, overlap=True)

    image.add_colorbar()  # unit is the nanomaggy, see http://www.sdss3.org/dr8/algorithms/magnitudes.php#nmgy
    image.colorbar.set_axis_label_font('serif', size=15)
    image.colorbar.set_axis_label_text(r'$\times$ 3.631 $\times$ 10$^{-6}$ Jy')
    image.colorbar.set_width(0.1)
    image.colorbar.set_font('serif', size=15)
    image.colorbar.set_pad(0.3)

    image.add_scalebar(radius / 4)
    image.scalebar.set_label('15"')
    image.scalebar.set_font('serif')
    image.scalebar.set_font_size(15)
    image.scalebar.set_color('black')

    image.axis_labels.set_xtext('Right ascension')
    image.axis_labels.set_ytext('Declination')
    image.set_axis_labels_family('serif')
    image.set_axis_labels_size(15)

    image.tick_labels.set_xformat('ddd.ddd')
    image.tick_labels.set_yformat('ddd.ddd')
    image.set_tick_labels_family('serif')
    image.set_tick_labels_size(15)


    image.add_label(0.83, 0.93, f'{levels[0] * 1e6:.0f} \N{GREEK SMALL LETTER MU}Jy beam\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}', family='serif', size=15, color='#e6194b', relative=True)

    save = f'/data5/sean/deep-fields/images/optical-{name}.png'
    image.save(save)
    print(f'Image saved! View it with gpicview {save}')
    # sys.exit()  # testing only


def main():
    blazars = '/data5/sean/deep-fields/catalogues/deep.fields.11.06.2019.cat.csv'
    df = pd.read_csv(blazars)

    for index, row in df.iterrows():  # get name, position, noise, optical image, and radio image
        optical = '/data5/sean/deep-fields/images/sdss-{}.fits'.format(row['Source name'])  # queried sdss manually

        if row['RA (J2000.0)'] > 156 and row['RA (J2000.0)'] < 168:  # lockman hole
            radio = '/data5/sean/deep-fields/catalogues/lockman.hole.11.06.2019.img.fits'
        elif row['RA (J2000.0)'] > 214 and row['RA (J2000.0)'] < 222:  # bootes
            radio = '/data5/sean/deep-fields/catalogues/bootes.11.06.2019.img.fits'
        else:  # elias n1, 237 < ra < 248
            radio = '/data5/sean/deep-fields/catalogues/elias.n1.11.06.2019.img.fits'

        do_plotting(name=row['Source name'], ra=row['RA (J2000.0)'], dec=row['Dec (J2000.0)'],
                    noise=row['Isl_rms'], optical=optical, radio=radio)


if __name__ == '__main__':
    main()
