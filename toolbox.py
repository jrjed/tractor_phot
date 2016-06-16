'''
A module containing some functions I found useful.
A few of which are used in the main fitting code
'''
from __future__ import division
import os
import sys

import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import montage_wrapper as montage
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord


def adxy(header, ra, dec):
    '''
    Converts RA and Dec to x and y pixel positions.
    '''
    ra = np.array(ra)
    dec = np.array(dec)
    x, y = WCS(header).all_world2pix(ra, dec, 0)
    return x, y


def angsep(ra0, dec0, ra, dec):
    '''
    Gives angular separation between the input position and the position after
    fitting in arcsec. Note: input should be in degrees.
    '''
    ra0 = np.array(ra0)
    dec0 = np.array(dec0)
    ra = np.array(ra)
    dec = np.array(dec)
    c1 = SkyCoord(ra0 * u.deg, dec0 * u.deg, frame='icrs')
    c2 = SkyCoord(ra * u.deg, dec * u.deg, frame='icrs')
    return c1.separation(c2).value * 3600.


def data2table(arrays=[], names=[], units=[], filename='d2table.fits'):
    '''
    Given a list of arrays, names and units, constructs table and writes
    to file. Units should be given as astropy.units.
    If extenstion of filename is .dat (.fits) it will be saved as an
    ascii (FITS) table.
    Returns the table that it saved to file.
    '''
    t = Table()
    if len(units) > 0:
        for info in zip(names, arrays, units):
            t.add_column(
                Column(name=info[0], data=info[1], unit=info[2]))
    else:
        for nm, arr in zip(names, arrays):
            t.add_column(Column(name=nm, data=arr))
    if os.path.splitext(filename)[1] == '.dat':
        t.write(filename, format='ipac')
    elif os.path.splitext(filename)[1] == '.fits':
        t.write(filename)
    else:
        print 'Unacceptable filename given. Table not being saved.'
    return t


def get_pixscale(header):
    pixscale = 0
    keys = header.keys()
    for key in keys:
        if 'PIXSCALE' in keys:
            pixscale = abs(header['PIXSCALE'])
        if 'PXSCAL1' in keys:
            pixscale = abs(header['PXSCAL1'])
        if 'CDELT1' in keys:
            pixscale = abs(header['CDELT1']) * 3600.  # in arcsec/pixel
        if 'CD1_1' in keys:
            pixscale = abs(header['CD1_1']) * 3600.  # in arcsec/pixel

    return pixscale


def nmatch(ra0, dec0, ra, dec, nneb):
    '''
    Returns indices of the nth nearest positional counterpart and the angular
    separation between them.
    '''
    ra0 = np.array(ra0)
    dec0 = np.array(dec0)
    ra = np.array(ra)
    dec = np.array(dec)
    reference_cat = SkyCoord(ra0 * u.deg, dec0 * u.deg, frame='icrs')
    other_cat = SkyCoord(ra * u.deg, dec * u.deg, frame='icrs')
    indices, angseps, physeps = reference_cat.match_to_catalog_sky(
        other_cat, nneb)
    return indices, angseps.to(u.arcsec).value


def matching(ra0, dec0, ra, dec):
    '''
    Returns indices of nearest positional counterpart and the angular
    separation between them.
    '''
    indices, angseps = nmatch(ra0, dec0, ra, dec, 1)
    return indices, angseps


def progress_bar(iterator, array, barLen=50):
    '''
    Displays progress bars for "for" loops.
    Input: i, array/list being looped over, optional: barLen = integer
    '''
    fraction_done = float(iterator + 1) / len(array)
    sys.stdout.write("\r")
    progress = ""
    for i in range(barLen):
        if i < int(barLen * fraction_done):
            progress += "#"
        else:
            progress += " "
    sys.stdout.write('Progress ')
    sys.stdout.write("[ %s ] %.0f %%" % (progress, fraction_done * 100))
    sys.stdout.flush()
    if iterator == len(array)-1:
        print " "


def read_image(image_path, extension=0):
    hdu = fits.open(image_path)
    image, header = hdu[extension].data, hdu[extension].header
    # To save some memory close the hdulist
    hdu.close()
    return image, header


def remove_existing(filename, save_directory, verbose=True):
    '''
    Overwrites file in the save directory that has a file name the same as
    the given filename.
    '''
    files = os.listdir(save_directory)
    if filename in files:
        os.remove(os.path.join(save_directory, filename))
    if verbose:
        print 'Removing ', str(filename)
    return


def res_sum(x0, y0, chi, half_fwhm_pix):
    '''
    Sums the residual flux in the beam around the position at which extraction
    took place.
    '''
    x_floor, y_floor = int(x0 + 0.5), int(y0 + 0.5)
    res = chi[y_floor - half_fwhm_pix: y_floor + half_fwhm_pix + 1,
              x_floor - half_fwhm_pix: x_floor + half_fwhm_pix + 1]
    left_over = res.sum()
    return left_over


def subimage(image_path, ra, dec, boxsize, showsub=False):
    '''
    Returns a subimage and subheader. Boxsize should be given in arcsec.
    This is to create a temporary subimage of a source for a quick view. To
    save the subimage, do it externally using
    astropy.io.fits.PrimaryHDU.writeto or montage_wrapper.mSubimage.
    '''
    _, header = read_image(image_path)
    x, y = adxy(header, ra, dec)
    # Center subimage on pixel source is located in (better centering)
    xpix, ypix = int(x + 0.5), int(y + 0.5)
    ra_center, dec_center = xyad(header, xpix, ypix)
    montage.commands.mSubimage(in_image=image_path,
                               out_image='temp.fits',
                               ra=ra_center,
                               dec=dec_center,
                               xsize=boxsize / 3600)
    subimage, subheader = read_image('temp.fits')
    remove_existing('temp.fits', os.getcwd(), verbose=False)

    if showsub:
        pdict = dict(interpolation='nearest',
                     origin='lower',
                     cmap='gray')
        plt.imshow(subimage, **pdict)
        plt.show()

    return subimage, subheader


def stamps2file(images=[], headers=[], names=[], file_path=None):
    '''
    Saves list of images, headers, and extension names to a FITS file
    '''
    primary = fits.PrimaryHDU()

    if len(headers) == 0:
        headers = [primary.header] * len(images)
    if len(names) == 0:
        names = ['none'] * len(images)

    hdus = [primary]
    for im, hdr, nm in zip(images, headers, names):
        hdu = fits.ImageHDU(data=im, header=hdr, name=nm)
        hdus.append(hdu)
    hdulist = fits.HDUList(hdus)
    hdulist.writeto(file_path, clobber=True)
    return hdulist


def xyad(header, x, y):
    '''
    Converts x and y pixel positions to RA and Dec.
    '''
    x = np.array(x)
    y = np.array(y)
    ra, dec = WCS(header).all_pix2world(x, y, 0)
    return ra, dec
