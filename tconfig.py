'''
Configuration file for tractor_phot
'''
import os

path2tractor_phot = os.path.dirname(os.path.abspath(__file__))

path2input = os.path.join(path2tractor_phot, 'input_files')

path2images = os.path.join(path2tractor_phot, 'images')

path2priors = os.path.join(path2tractor_phot, 'priors_catalogs')

path2output = os.path.join(path2tractor_phot, 'output')

path2psfs = os.path.join(path2tractor_phot, 'psfs')
