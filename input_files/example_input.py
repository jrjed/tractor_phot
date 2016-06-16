'''
TEST INPUT FILE (IRAC CHANNEL 1)
'''
# Path Information
image_filename= 'example_image.fits' # image filename
image_ext = 0 # The extension of the image in the fits file

# Prior Information
priors_filename = 'example_priors.fits' # priors fits table filename

# Image Information
sky_noise =  0.002489 # std from several empty apertures on the map
sky_level =  0.09841  # mean from several empty apertures on the map

# PSF Info
FWHM = 1.95 # in arcsec from IRAC Instrument Handbook (warm value)
psf_filename = 'example_psf.fits' # optional set to None if not used

# Fitting Logistics
iterations = 'all' # set to 'all' to run for all sources in priors catalog
boxsize = 6 # length of side of subimages in arcsec

# Output
output_directory = 'example' # The directory where output should be saved.
save_params = True
save_stamps = True

'''
PERSONAL NOTES:
    This is an example/test input file.

'''
