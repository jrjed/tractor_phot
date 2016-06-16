'''
Plots a quick comparison between the Tractor fluxes and catalog input.

A more robust comparison than this may be required for your results.

This is run from your shell command line.

Example use:
    python compare.py irac1 b1_kpriors.fits point_output.fits

To use Seaborn:
    python compare.py -pretty irac1 b1_kpriors.fits point_output.fits

This code assumes that priors with a flux of 1.0 do not have a catalog
counterpart to compare too.
'''
import os
import sys
import tconfig
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from astropy.table import Table

def filter_table(in_table, out_table):
    '''
    Filter out invalid values from the FITS tables
    '''
    filt = np.where((in_table['flux'] != 1.0) & \
                    (out_table['flux'] > 0) & \
                    (out_table['mdiff'] != -99))[0]
    #nans will be automatically excluded from mean and std calculations
    return in_table[filt], out_table[filt]

def flux2ab(flux, conversion=1):
    '''
    Converts flux in map units to AB magnitudes
    '''
    flux *= conversion
    return -2.5 * np.log10(flux.data/(3631.e6))

def comp_stats(mdiff):
    '''
    Compute stats to be displayed on the plot.
    '''
    stats = []
    dround = 4

    mu = np.nanmean(mdiff)
    sigma = np.nanstd(mdiff)

    stats.append(r'$\Delta$m = {}'.format(round(mu, dround)))
    stats.append(r'$\sigma$ = {}'.format(round(sigma, dround)))
    return '\n'.join(stats)

def plot_results(mdiff, mcat, stats, channel='', blend=[]):
    '''
    Plot the offset vs. catalog magnitude
    '''
    mu = np.nanmean(mdiff)
    sigma = np.nanstd(mdiff)

    fig, ax = plt.subplots(1,1)
    text0 = AnchoredText(stats, loc=2)
    ax.set_title(channel.upper())
    ax.scatter(mcat, mdiff)
    #ax.hexbin(mcat, mdiff, cmap='inferno')
    ax.axhline(0, color='k')
    ax.axhline(mu, linestyle='--', color='k', label=r'$\Delta$m')
    ax.axhline(mu + sigma, linestyle='--', color='r', label=r'$\sigma$')
    ax.axhline(mu - sigma, linestyle='--', color='r')
    ax.set_xlabel('$m_{Catalog}$', fontsize=14)
    ax.set_ylabel('$m_{Tractor}$ - $m_{Catalog}$', fontsize=14)
    ax.add_artist(text0)

    if len(blend) > 0:
        ax.scatter(mcat[blend], mdiff[blend], color='c', label='Blended/Crowded')

    ax.legend(loc='best')
    plt.show()
    #plt.savefig('/Users/Jesse/projects/eros/writeups/' +channel+'_comp.png')

def main():
    #Take in arguments
    if len(sys.argv) > 4:
        if sys.argv[1].lower() == '-pretty':
            import seaborn
            channel = sys.argv[2].lower()
            infile = os.path.join(tconfig.path2priors, sys.argv[3])
            outfile = os.path.join(tconfig.path2output, channel, sys.argv[4])
        else:
            print 'Invalid fourth argument!'
    else:
        channel = sys.argv[1].lower()
        infile = os.path.join(tconfig.path2priors, sys.argv[2])
        outfile = os.path.join(tconfig.path2output, channel, sys.argv[3])

    #Read in the FITS tables
    tout = Table.read(outfile)
    tin = Table.read(infile)[:len(tout)]
    tin, tout = filter_table(tin, tout)
    print '{} sources included in comparison'.format(len(tout))

    #Figure out conversion factor to uJy
    if channel in ['irac1', 'irac2', 'irac3', 'irac4']:
        flux_conv = 8.461594994075238 #assumes 0.6'' pixel size
        sr = 3
    elif channel in ['mips1']:
        flux_conv = 33.84637997630095 #assumes 1.2'' pixel size
        sr = 6
    else:
        print 'Band not recognized. Leaving in map units * (pix**2)'
        flux_conv = 1.0

    if channel == 'irac1':
        flux_conv *= 1.02 # additional calibration correction factor

    #Convert fluxes to magnitudes
    mc = flux2ab(tin['flux'], flux_conv) #should be in map units e.g. (MJy/sr) * pixel^2
    mt = flux2ab(tout['flux'], flux_conv)
    mdiff = mt - mc

    #Calculate stats to put on the plot
    stats = comp_stats(mdiff)
    blend = []

    #Testing overplot blends. Uncomment this block to overplot blended sources
    import myth.cerca as cc
    infile = os.path.join(tconfig.path2priors, sys.argv[3])
    t = Table.read(infile)
    ra = t['ra']
    dec = t['dec']
    ra0 = tout['ra']
    dec0 = tout['dec']
    idx, sep = cc.nmatch(ra0, dec0, ra, dec, 2)
    blend = np.where(sep <= sr)[0]
    print '{} sources are blended'.format(len(blend))

    plot_results(mdiff, mc, stats, channel, blend)

if __name__ == '__main__':
    main()

