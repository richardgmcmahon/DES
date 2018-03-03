"""

 read a DESDM SExtractor psfcat file and do some analysis

 need to read the PSFEx config file to get some the default parameters

"""
# make future proof for Python 3
from __future__ import print_function, division

import os
import sys
import time
from time import strftime, gmtime

import getpass
import inspect
import logging
import traceback

# read in various libraries and print the version numbers
import matplotlib as mpl
print('matplotlib: ', mpl.__version__)
import matplotlib.pyplot as plt

import numpy as np
print('numpy: ', np.__version__)

import scipy
print('scipy: ', scipy.__version__)
from scipy.ndimage import interpolation

import astropy
print('astropy: ', astropy.__version__)
from astropy.table import Table, join
from astropy.io import fits
from astropy.stats import median_absolute_deviation

sys.path.append('/home/rgm/soft/python/lib/')
import psfex
try:
    print('psfex: ', psfex.__version__)
except:
    pass

from librgm.plotid import plotid

import psfex_util as psfex_util

from radial_profile import *
print('inspect.getfile(radial_profile_greenfield): ',
    inspect.getfile(radial_profile_greenfield))


__version__ = 0.1

def psf_psfcat_plots(infile_psfcat=None,
                     infile_psfex_starlist=None,
                     RELEASE='Y1A1', TILE=None,
                     WAVEBAND='i', pause=False):
    """

    """

    if infile_psfcat is None:
        datapath='/data/desardata/'+ RELEASE + '/' + TILE + '/'

        datapath='/data/desardata/PSFEX/'+ RELEASE + '/' + TILE + '/R1/'

        filename_image= TILE+ '_' + WAVEBAND + '.fits.fz'
        infile_image= datapath + filename_image

        filename_psf= TILE + '_' + WAVEBAND + '_psfcat.psf'
        infile_psf = datapath + filename_psf

        #if RELEASE == 'SVA1':
        #if RELEASE == 'Y1A1':

        infile_psfcat= datapath + TILE + '_' + WAVEBAND + '_psfcat.fits.fz'

        infile_psfcat= datapath + TILE + '_' + WAVEBAND + '_psfcat.fits'

        infile_psfex_starlist = datapath + TILE + '_' + WAVEBAND + '_psfcat-starlist.fits'

    print('Reading in: ', infile_psfex_starlist)
    starlist = rd_psfex_starlist(infile=infile_psfex_starlist, pause=pause)
    print('Read in: ', infile_psfex_starlist)

    starlist.info('stats')
    trace = traceback.extract_stack()
    print('traceback: ',trace[0][0], '; ',trace[0][1])
    if pause:
        key=raw_input("Enter any key to continue: ")

    print('Reading in: ', infile_psfcat)
    psfcat = rd_psfcat(infile=infile_psfcat, plots=False, pause=pause)
    print('Read in: ', infile_psfcat)

    psfcat.info('stats')
    trace = traceback.extract_stack()
    print('traceback: ',trace[0][0], '; ',trace[0][1])
    if pause: key=raw_input("Enter any key to continue: ")

    # plot some PSFEx diagnostics
    psfex_util.psfcat_plots(data=psfcat, infile=infile_psfcat,
        label= RELEASE,
        overplot=False,
        RANGE_FLUX_RADIUS=(-0.5, 6.5),
        RANGE_FLUX_APER=(10.0, 1E7),
        RANGE_ELONGATION=(0.5, 4.0))

    # now join

    # now match with the psfcat file with starlist
    # could write it out as a FITS file
    psfcat_starlist = match_starlist(psfcat=psfcat, starlist=starlist)

    psfex_util.psfcat_plots(data=psfcat_starlist, infile=infile_psfcat,
        label= RELEASE + '_' + 'starlist',
        overplot=True,
        RANGE_FLUX_RADIUS=(-0.5, 6.5),
        RANGE_FLUX_APER=(10.0, 1E7),
        RANGE_ELONGATION=(0.5, 4.0))

    #sys.exit()



    # plot zoom in plots if you need them
    plot_zoom = True
    if plot_zoom:

        plt.close()
        markersize = 4
        # zoom in
        psfex_util.psfcat_plots(data=psfcat, infile=infile_psfcat,
            label= RELEASE + '_' + 'zoom',
            overplot=False,
            markersize=markersize,
            RANGE_FLUX_RADIUS=(1.0, 4.0),
            RANGE_FLUX_APER=(100.0, 1E6))

        psfex_util.psfcat_plots(data=psfcat_starlist, infile=infile_psfcat,
            label= RELEASE + '_' + 'starlist_zoom',
            overplot=True,
            markersize=markersize,
            RANGE_FLUX_RADIUS=(1.0, 4.0),
            RANGE_FLUX_APER=(100.0, 1E6))



if __name__ == '__main__':
    """

    """

    import argparse
    import ConfigParser

    t0 = time.time()
    datestamp = time.strftime("%Y%m%d", gmtime())

    release_default = 'Y1A1'
    tile_default = 'DES2327-5248'
    file_default = None

    pause=False

    RELEASE = 'SVA1'
    RELEASE = 'Y1A1'

    WAVEBAND='g'
    WAVEBAND='i'

    WAVEBANDS = ['g','r','i','z','Y']

    TILE='DES1000+0209'

    #datapath = datapath_root + '/' + RELEASE + '/COSMOS/' + TILE + '/'

    TILE = 'DES0449-4706'
    TILE = 'DES0449-4748'


    TILE = 'DES2327-5248'

    TILE = 'DES0406-5414'

    description = 'DES psfcat analysis'
    epilog = ''
    parser =  argparse.ArgumentParser(
        description=description, epilog=epilog,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--release", default=release_default,
        dest='release', help="DESDM Release as a string list e.g. Y1A1 ")

    parser.add_argument("--tile", default=tile_default,
        dest='tile', help="DESDM Tile e.g. DES0000+0000 ")

    parser.add_argument("--infile_psfcat", default=None,
        help="psfex input psfcat SExtractor vignette catalogue file")

    parser.add_argument("--infile_psfex_starlist", default=None,
        help="psfex output starlist file")

    parser.add_argument("--format", default='coadd',
        dest='format', help="coadd or SE")

    parser.add_argument("--version", action="store_true", dest="version",
                  default="", help="print version number and  exit")

    parser.add_argument("--pause", action="store_true",
        dest='pause', help="turn of pausing option")

    print('Number of arguments:', len(sys.argv), 'arguments: ',sys.argv[0])

    args = parser.parse_args()

    if args.version:
        print('version: ',__version__)
        sys.exit(0)

    pause = False
    if args.pause:
        pause = True

    print('HOME: ', os.environ['HOME'])
    print('USERNAME: ',  getpass.getuser())
    config_path = os.environ['HOME'] + '/config/des/'
    # look in current directory first
    # check if it exists as give a warning or could use '/home/rgm/config/des/'

    format = args.format

    infile_psfcat = args.infile_psfcat
    if infile_psfcat is not None:
        print('Procssing file:', infile_psfcat)

    infile_psfex_starlist = args.infile_psfex_starlist
    if infile_psfex_starlist is not None:
        print('Procssing file:', infile_psfex_starlist)


    TILE = args.tile
    print('Processing tile: ', TILE)
    RELEASE = args.release
    print('Processing release: ', RELEASE)

    Config = ConfigParser.RawConfigParser()

    # read config file
    Config.read('psfcat_analysis.cfg')
    print('Config.sections(): ', Config.sections())

    datapath = '/data/desardata/SVA1/COSMOS/' + TILE + '/'

    datapath_root = Config.get('local','datapath_root')
    print('datapath_root: ', datapath_root)


    apertures = set_des_apertures(unit='pixels')

    wavebands=['g','r','i','z','Y']

    t0=time.time()


    if str.lower(format) == 'coadd':

        for WAVEBAND in WAVEBANDS:

            psf_psfcat_plots(RELEASE='Y1A1', TILE=TILE, WAVEBAND=WAVEBAND)

            fignums = plt.get_fignums()
            print(fignums)
            print('plt.get_fignums(): ', plt.get_fignums())

            for ifig in fignums:
                plt.close(ifig)

            print('plt.get_fignums(): ', plt.get_fignums())

            key=raw_input("Enter any key to continue: ")

    if str.upper(format) == 'SE':

            plt.close()
            psf_psfcat_plots(
                infile_psfcat=infile_psfcat,
                infile_psfex_starlist=infile_psfex_starlist)

            fignums = plt.get_fignums()

            print(fignums)
            print('plt.get_fignums(): ', plt.get_fignums())

            for ifig in fignums:
                plt.close(ifig)

            print('plt.get_fignums(): ', plt.get_fignums())

            print('Reading starlist:', infile_psfex_starlist)
            psfex_starlist = Table.read(infile_psfex_starlist,  hdu=2)
            psfex_starlist.info('stats')

            xdata = psfex_starlist['X_IMAGE']
            ydata = psfex_starlist['Y_IMAGE']
            print('x range:', np.min(xdata), np.max(xdata))
            print('x range:', np.min(ydata), np.max(ydata))
            ndata = len(xdata)
            plt.plot(xdata, ydata, '.', label=str(ndata))
            plt.title(infile_psfex_starlist, fontsize='small')
            plt.legend()

            plt.xlim(0.0, 2048.0)
            plt.ylim(0.0, 4096.0)
            plt.xlabel('X_IMAGE')
            plt.ylabel('Y_IMAGE')
            plt.axes().set_aspect('equal')
            plt.show()

            key=raw_input("Enter any key to continue: ")
