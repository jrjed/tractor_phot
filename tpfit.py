#!/usr/bin/python
'''
This code fits a given PSF model at the given positions using the given priors.
'''
from __future__ import division
import os
import sys
import time
import subprocess

import numpy as np
from astropy.table import Table
from tractor.image import Image
from tractor.wcs import PixPos, NullWCS
from tractor.psf import PixelizedPSF, NCircularGaussianPSF
from tractor.sky import ConstantSky
from tractor.brightness import Flux, NullPhotoCal
from tractor.engine import Tractor
from tractor.pointsource import PointSource

import tconfig
import toolbox as tb


class ScImage(object):
    '''
    Contains image properties and methods related to science images
    '''

    def __init__(self, image_file, extension=0, noise=None, level=None,
                 fwhm=None, psf_file=None, boxsize=None):
        self.image_file = image_file
        self.image_path = os.path.join(tconfig.path2images, image_file)
        self.extension = extension
        self.noise = noise
        self.level = level
        self.fwhm = fwhm
        self.psf_file = psf_file
        self.image, self.header = tb.read_image(self.image_path,
                                                self.extension)
        self.pixscale = tb.get_pixscale(self.header)

        if boxsize is None:
            self.boxsize = np.ceil(3 * self.fwhm)
        else:
            self.boxsize = boxsize  # size of subimage in arcsec to be cutout

    def subimage(self, ra, dec):
        '''
        Returns a subimage and subheader. Boxsize should be given in arcsec
        Uses Montage + montage_wrapper
        '''

        subimage, subheader = tb.subimage(self.image_path,
                                          ra,
                                          dec,
                                          self.boxsize)
        return subimage, subheader

    def psf_model(self):
        '''
        Determines the PSF model to use
        '''
        if self.psf_file is not None:
            psf_image, _ = tb.read_image(
                os.path.join(tconfig.path2psfs, self.psf_file))
            psfmod = PixelizedPSF(psf_image)
        else:
            psf_sigma = self.fwhm / np.sqrt(8 * np.log(2)) / self.pixscale
            psfmod = NCircularGaussianPSF([psf_sigma], [1.0])
        return psfmod

    def tractor_image(self, subimage):
        '''
        Creates the Tractor image
        '''
        timage = [Image(data=subimage,
                        invvar=np.ones_like(subimage) / self.noise**2,
                        psf=self.psf_model(),
                        wcs=NullWCS(),
                        photocal=NullPhotoCal(),
                        sky=ConstantSky(self.level))]
        return timage

    def inbox(self, subheader, priors):
        '''
        Returns indices of sources that lie inside the given subimage
        '''
        x, y = tb.adxy(subheader, priors['ra'], priors['dec'])

        # Only need xc as xc = yc for a square subimage
        xc = int(subheader['NAXIS1'] / 2 - 0.5)

        halfbox = int(np.ceil(self.boxsize / 2 / self.pixscale))
        extra = self.fwhm / self.pixscale
        search_length = halfbox + 0.5 + extra
        boxlim = [xc - search_length, xc + search_length]
        inbox = np.where((x >= boxlim[0]) & (x <= boxlim[1]) &
                         (y >= boxlim[0]) & (y <= boxlim[1]))[0]

        return priors[inbox]


def get_priors(priors_filename):
    '''
    Reads the priors catalog into a table
    '''
    priors_path = os.path.join(tconfig.path2priors, priors_filename)
    table = Table.read(priors_path)

    if 'flux' not in table.keys():
        table['flux'] = np.ones(len(table))

    if 'id' not in table.keys():
        table['id'] = ['source'+str(i+1) for i in xrange(len(table))]

    return table


def pointsource_fit(xsub, ysub, fluxsub, TractorImage):
    '''
    Returns the optimized Tractor object and variances of each source fit
    '''
    # Generate Tractor source object
    TractorSources = [PointSource(PixPos(x, y),
                                  Flux(flux))
                      for x, y, flux in zip(xsub, ysub, fluxsub)]

    # Freeze position
    for src in TractorSources:
        src.pos.freezeParam('x')
        src.pos.freezeParam('y')

    # Construct Final Tractor Object for Fitting
    tractor = Tractor(TractorImage, TractorSources)

    # Freeze image parameters
    tractor.freezeParam('images')

    # Initial Source Model and Initial Residual left behind in map units
    # mod0 = tractor.getModelImage(0)
    # res0 = tractor.getChiImage(0) * sky_noise

    # Fit the model
    for k in xrange(100):
        dlnp, X, alpha, var = tractor.optimize(variance=True)
        # print 'dlnp',dlnp
        # print 'var', var
        if dlnp < 1e-3:
            break
    return tractor, var


def main():

    # Timing
    start = time.time()

    # Take in parameters from input file
    input_file = str(os.path.join(tconfig.path2input, sys.argv[1]))
    subprocess.call(['cp', input_file, 'tempinput.py'])
    import tempinput
    subprocess.call(['rm', 'tempinput.py'])
    subprocess.call(['rm', 'tempinput.pyc'])

    image_filename = tempinput.image_filename
    image_ext = tempinput.image_ext
    output_directory = tempinput.output_directory
    priors_filename = tempinput.priors_filename
    sky_noise = tempinput.sky_noise
    sky_level = tempinput.sky_level
    FWHM = tempinput.FWHM
    psf_filename = tempinput.psf_filename
    iterations = tempinput.iterations
    boxsize = tempinput.boxsize
    save_params = tempinput.save_params
    save_stamps = tempinput.save_stamps

    # Read in the Science Image with parameters set in input file
    science = ScImage(image_file=image_filename,
                      extension=image_ext,
                      noise=sky_noise,
                      level=sky_level,
                      fwhm=FWHM,
                      psf_file=psf_filename,
                      boxsize=boxsize)

    # Read in the source priors
    priors = get_priors(priors_filename)

    # Determine exact x,y position in full science image
    xtot, ytot = tb.adxy(science.header, priors['ra'], priors['dec'])

    # Create empty arrays to fill with values
    if iterations == 'all':
        iterations = len(priors)
    elif iterations > len(priors):
        print 'ERROR: ITERATIONS IS LARGER THEN NUMBER OF SOURCES IN TABLE'

    fluxout = np.empty(iterations)
    rms = np.empty(iterations)
    mdiff = np.empty(iterations)
    flux_err = np.empty(iterations)

    for i, prior in enumerate(priors[:iterations]):

        # Create subimage of the source
        subimage, subheader = science.subimage(prior['ra'],
                                               prior['dec'])

        # Generate Tractor Image
        TractorImage = science.tractor_image(subimage)

        # Find all sources inside and just outside the subimage
        subpriors = science.inbox(subheader, priors)
        xsub, ysub = tb.adxy(subheader, subpriors['ra'], subpriors['dec'])
        fluxsub = subpriors['flux']
        central_source = np.where(subpriors['id'] == prior['id'])[0]

        # Fit the sources in the subimage
        tractor, var = pointsource_fit(xsub,
                                       ysub,
                                       fluxsub,
                                       TractorImage)

        # Optimized model and Residaul left behind in map units
        mod = tractor.getModelImage(0)
        res = tractor.getChiImage(0) * science.noise

        # Fill output arrays
        fluxout[i] = tractor.getParams()[central_source]
        rms[i] = res.std()
        flux_err[i] = np.sqrt(var[central_source])

        if prior['flux'] == 1.0:
            mdiff[i] = np.nan
        else:
            mdiff[i] = -2.5 * np.log10(fluxout[i] / prior['flux'])

        if save_stamps:
            stamps = [subimage, mod, res]
            headers = [subheader] * 3
            snames = ['science', 'model', 'residual']
            stamp_path = os.path.join(tconfig.path2output,
                                      output_directory,
                                      prior['id'] + '.fits')

            tb.stamps2file(stamps, headers, snames, stamp_path)

        tb.progress_bar(i, xrange(iterations))

        if i == (iterations - 1):
            break

    if save_params:
        priors = priors[:iterations]
        priors['flux'] = fluxout
        priors['residual_rms'] = rms
        priors['mdiff'] = mdiff
        priors['flux_err'] = flux_err
        priors = priors['id',
                        'ra',
                        'dec',
                        'flux',
                        'flux_err',
                        'residual_rms',
                        'mdiff']  # reorder the table for neatness

        # Save/overwrite output table
        print 'Saving/overwriting tracphot_output.fits'
        tb.remove_existing(os.path.join(tconfig.path2output,
                                        output_directory,
                                        'tracphot_output.fits'),
                           verbose=False)

        priors.write(os.path.join(tconfig.path2output,
                                  output_directory,
                                  'tracphot_output.fits'))

    end = time.time()
    print '{} minutes to complete'.format(round((end-start)/60, 1))

    if sys.argv[1] == 'example_input.py':
        print 'TEST SUCESSFUL!'

if __name__ == '__main__':
    main()
