#!/usr/bin/env python3

''' Get images of the blazars in the Deep Fields data. '''

import matplotlib as mpl
mpl.use('Agg')

import argparse
import aplpy
import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy import coordinates
from astropy import units as u

__author__ = 'Sean Mooney'
__email__ = 'sean.mooney@ucdconnect.ie'
__date__ = '11 June 2019'

def make_cut_out_image(sources, radius=1 / 60, cmap='viridis', vmin=0, output='', contours=[5, 10, 20, 40]):
    ''' Make a cut-out image of a given source. '''

    df = pd.read_csv(sources)

    for source_name, ra, dec, peak_flux, noise in zip(df['Source name'], df['RA'], df['DEC'], df['Peak_flux'], df['Isl_rms']):
        if ra > 156 and ra < 168:
            field_file = '/data5/sean/deep-fields/catalogues/lockman.hole.11.06.2019.img.fits'
        if ra > 214 and ra < 222:
            field_file = '/data5/sean/deep-fields/catalogues/bootes.11.06.2019.img.fits'
        if ra > 237 and ra < 248:
            field_file = '/data5/sean/deep-fields/catalogues/elias.n1.11.06.2019.img.fits'

        image = aplpy.FITSFigure(field_file)  # BUG running this on my laptop gives a MemoryError

        # image.show_regions('/data5/sean/deep-fields/ellipses.reg')
        image.recenter(ra, dec, radius=radius)  # BUG see https://github.com/aplpy/aplpy/issues/423#issuecomment-478535924
        image.show_colorscale(cmap=cmap, vmin=vmin, vmax=peak_flux, stretch='arcsinh')
        image.show_contour(levels=[contour * noise for contour in contours], colors=['#f8b1af', '#f4827f', '#f0544f', '#af3e3a'], overlap=True, alpha=0.67)

        image.add_colorbar()
        image.colorbar.set_axis_label_font('serif', size=15)
        image.colorbar.set_axis_label_text(r'Jy beam$^{-1}$')
        image.colorbar.set_width(0.1)
        image.colorbar.set_font('serif', size=15)
        image.colorbar.set_pad(0.3)

        image.add_scalebar(radius / 4)
        image.scalebar.set_label('15"')
        image.scalebar.set_font('serif')
        image.scalebar.set_font_size(15)
        image.scalebar.set_color('white')

        image.add_beam()
        image.beam.set_linestyle('dashed')
        image.beam.set_edgecolor('white')
        image.beam.set_facecolor('none')

        image.axis_labels.set_xtext('Right ascension')
        image.axis_labels.set_ytext('Declination')
        image.set_axis_labels_family('serif')
        image.set_axis_labels_size(15)

        image.tick_labels.set_xformat('ddd.ddd')
        image.tick_labels.set_yformat('ddd.ddd')
        image.set_tick_labels_family('serif')
        image.set_tick_labels_size(15)

        image.add_label(0.83, 0.93, f'\N{GREEK SMALL LETTER SIGMA} = {noise * 1e6:.0f} \N{GREEK SMALL LETTER MU}Jy beam\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}', family='serif', size=15, color='white', relative=True)

        image.save('/data5/sean/deep-fields/images/{}.png'.format(source_name))
        print('Done! View it: gpicview /data5/sean/deep-fields/images/{}.png'.format(source_name))
        # sys.exit()
        # TODO get colorbar ticks on the inside


def main():
    formatter_class = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=formatter_class)

    parser.add_argument('-s', '--sources', required=False, type=str, help='CSV of the blazars', default='/data5/sean/deep-fields/catalogues/deep.fields.11.06.2019.cat.csv')

    args = parser.parse_args()
    sources = args.sources
    make_cut_out_image(sources)


if __name__ == '__main__':
    main()
