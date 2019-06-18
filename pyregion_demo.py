import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import WCSAxes
import pyregion


# Create figure
fig = plt.figure(figsize=(8, 4))

# Parse WCS information
hdu = fits.open('/mnt/closet/ldr2-blazars/deep-fields/bootes-image.fits')[0]
wcs = WCS(hdu.header)

# Create axes
import numpy as np
asdf = np.array([[0.1], [0.1], [0.4], [0.8]])
print(asdf.shape)
asdf = np.array([[[[0.1, 0.1, 0.4, 0.8]]]])
asdf=np.squeeze(asdf)
print(asdf.shape)
ax = WCSAxes(fig, asdf, wcs=wcs)
fig.add_axes(ax)

# Hide labels on y axis
ax2.coords[1].set_ticklabel_position('')

ax.set_xlim(300, 1300)
ax.set_ylim(300, 1300)
ax.set_aspect(1)

r = pyregion.open('/mnt/closet/ldr2-blazars/deep-fields/bootes-image.fits/ellipses.reg').as_imagecoord(header)

patch_list, text_list = r.get_mpl_patches_texts()

for p in patch_list:
    ax.add_patch(p)

ax.add_artist(atext)

plt.draw()
plt.show()
