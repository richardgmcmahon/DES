"""

 read a DESDM psfcat file and do some analysis

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

from psfex_util import *

from radial_profile import *
print('inspect.getfile(radial_profile_greenfield): ', 
    inspect.getfile(radial_profile_greenfield))


__version__ = 0.1

def psf_psfcat_plots(RELEASE='Y1A1', TILE=None, WAVEBAND='i', pause=False):
    """

    """

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
    if pause: key=raw_input("Enter any key to continue: ")

    print('Reading in: ', infile_psfcat)
    psfcat = rd_psfcat(infile=infile_psfcat, plots=False, pause=pause)
    print('Read in: ', infile_psfcat)

    psfcat.info('stats')  
    trace = traceback.extract_stack()
    print('traceback: ',trace[0][0], '; ',trace[0][1])
    if pause: key=raw_input("Enter any key to continue: ")

    # plot some PSFEx diagnostics 
    psfcat_plots(data=psfcat, infile=infile_psfcat,
        label= RELEASE,
        overplot=False, 
        RANGE_FLUX_RADIUS=(-0.5, 6.5),
        RANGE_FLUX_APER=(10.0, 1E7),
        RANGE_ELONGATION=(0.5, 4.0))

    # now join 

    # now match with the psfcat file with starlist 
    # could write it out as a FITS file
    psfcat_starlist = match_starlist(psfcat=psfcat, starlist=starlist)

    psfcat_plots(data=psfcat_starlist, infile=infile_psfcat, 
        label= RELEASE + '_' + 'starlist', 
        overplot=True,
        RANGE_FLUX_RADIUS=(-0.5, 6.5),
        RANGE_FLUX_APER=(10.0, 1E7),
        RANGE_ELONGATION=(0.5, 4.0))
 
    #sys.exit()

    # plot zoom in plots if you need them
    plot_zoom = False
    if plot_zoom:

        # zoom in   
        psfcat_plots(data=psfcat, infile=infile_psfcat,
            label= RELEASE + '_' + 'zoom',
            overplot=False, 
            RANGE_FLUX_RADIUS=(-0.5, 6.5),
            RANGE_FLUX_APER=(10.0, 1E7))

        psfcat_plots(data=psfcat_starlist, infile=infile_psfcat, 
            label= RELEASE + '_' + 'starlist_zoom', 
            overplot=True, 
            RANGE_FLUX_RADIUS=(-0.5, 6.5),
            RANGE_FLUX_APER=(10.0, 1E7))





if __name__ == '__main__':
    """

    """

    from argparse import ArgumentParser    
    import ConfigParser

    
    t0 = time.time()
    datestamp = time.strftime("%Y%m%d", gmtime())


    release_default = 'Y1A1'
    tile_default = 'DES2327-5248'

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



    parser = ArgumentParser(description=
      'DESDM psfcat analysis')


    parser.add_argument("--release", default=release_default,
        dest='release', help="DESDM Release as a string list e.g. Y1A1 ")


    parser.add_argument("--tile", default=tile_default,
        dest='tile', help="DESDM Tile e.g. DES0000+0000 ")

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


    for WAVEBAND in WAVEBANDS:
    
        psf_psfcat_plots(RELEASE='Y1A1', TILE=TILE, WAVEBAND=WAVEBAND)

        fignums = plt.get_fignums()
        print(fignums)
        print('plt.get_fignums(): ', plt.get_fignums())
    
        for ifig in fignums:
            plt.close(ifig)

        print('plt.get_fignums(): ', plt.get_fignums())

        key=raw_input("Enter any key to continue: ")

