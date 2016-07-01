
# tractor_phot #
tractor_phot is an implementation of [the Tractor](https://github.com/dstndstn/tractor) for PSF fitting unresolved sources. You can clone this repo and use as is and/or mangle it to your liking.

## Getting Started ##
tractor_phot is simply a directory from which the main fitting script runs. The code fits a given PSF model to the given coordinates.  Most of the installation lies in aquiring the required libraries.

Simply clone the repo,

```
$ git clone https://github.com/jrjed/tractor_phot
```
and add the tractor_phot directory to your PYTHONPATH.

### Requirements ###
The code requires *Python 2.7+*, as well as a few additional packages.

***[the Tractor](https://github.com/dstndstn/tractor)***

***[astropy](https://github.com/astropy/astropy)***

***[Montage](http://montage.ipac.caltech.edu/)***

***[montage-wrapper](http://www.astropy.org/montage-wrapper/)***

## Use ##

###1.) Input Files ###
An input file provides the logistical information that the code needs to run. 
An example input file ("example_input.py") can be found in the "input_files/" directory. It is here that the user will provide the following info:

***Science Image Properties***  
image_filename: filename of science image located in the "images/" directory  
image_ext: FITS image extension (defaults to 0)  
sky_noise: the noise density (e.g. per-pixel noise) of the science image in the map units  
sky_level: sky background of the science image in the map units  

***Priors***  
priors_filename: filename of priors FITS table in the "priors_catalogs/" directory  

***PSF Info***  
FWHM: instrument FWHM in arcseconds  
psf_filename (optional): filename of FITS image in "psfs/" to be used as PSF model. The PSF image should be rotated and interpolated to the header specifications of the science image. If not given (i.e. set to None), a circular Gaussian with the given FWHM will be used as the PSF model. See the example PSF in the "psfs/" directory.

***Fitting Logistics***  
iterations: number of sources to be fit e.g. if you just want to fit the first 100 sources (that is, if your primary sources are the first 100 sources in your priors catalog and the rest are peripheral sources) set to 100. Set to 'all' for all prior sources to be fit and included in the output.  
boxsize: length of side of subimages in arcsec. If set to None, the default is 3 X FWHM.  

***Ouput***  
output_directory: subdirectory of "output/" directory where the output should be saved  
save_params: save output to a fits table (i.e. tracphot_output.fits). Set True or False  
save_stamps: save the subimage, model, and residual in a single FITS image file. Set True or False  

###2.) Priors Catalog###

The priors catalog should be a FITS table containing prior knowledge of the sources to be fit. The code does not have any source detection feature and so all sources that may end up inside or just outside a subimage centered on the primary source should be included for proper cleaning of the subimage and/or deblending.

The script requires prior knowledge of each source's position. Optionally, one can also give prior fluxes, and source identifications. Source flux priors can be obtained from other catalogs and are useful for comparison to the Tractor extracted flux. Sources for which the user does not have flux priors should be initialized to unity. If no flux column is given in the priors table, all flux priors will automatically be initialized to unity. Additionally, if no source identification is given, sources are labeled as "source+n" in the output FITS table. 

The column headers should be named: "id", "ra", "dec", "flux". Of course, leave out "id" and "flux" columns as needed. The "data2table" function in "toolbox" makes it easy to save the desired data in the proper format. See "example_priors.fits" in the "priors_catalog/" directory for an example.

###3.) Running the code###

You'll want to put your input file in the , "input_files/" directory, priors FITS catalog in the "priors_catalogs/" directory, science image in the "images/" directory, and the (optional) PSF image in the "psfs/" directory. Additionally, add a subdirectory to the "output/" directory e.g. "irac1". This is where the output will be saved as indicated by your input file. 

The script is run from the command line as follows:
```
$ cd tractor_phot
$ python tpfit.py example_input.py
```

##Tests##
You can run a test to make sure things are working properly:
```
$ cd tractor_phot
$ python test.py
```
If working correctly, "TEST SUCESSFUL!" should print when done running (it should take a couple moments). Additionally, this should write a output FITS table and FITS images to the "output/example/" directory.
