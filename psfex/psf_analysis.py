"""
 read a DESDM psfcat starlist file and do some analysis of the psf and starlist
 vignettes

"""
# make future proof for Python 3
from __future__ import print_function, division

import os
import sys
import time
import inspect

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

import skimage
print('scikit-image: ', skimage.__version__)


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


if __name__ == '__main__':
    """

    """
    apertures = set_des_apertures(unit='pixels')

    Tile='DES0449-4748'

    wavebands=['g','r','i','z','Y']


    t0=time.time()

    pause=False

    RELEASE = 'SVA1'
    RELEASE = 'Y1A1'

    WAVEBAND='g'
    WAVEBAND='i'




    TILE='DES1000+0209'
    datapath='/data/desardata/SVA1/COSMOS/' + TILE + '/'

    TILE='DES0449-4706'
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

    # select some star to analyze; e.g. 10x10 brightest; a random set;
    # or slice in mag as a function of x, y location
    psfcat_starlist = match_starlist(psfcat=psfcat, starlist=starlist)


    psfcat_starlist.info('stats')
    trace = traceback.extract_stack()
    print('traceback: ',trace[0][0], '; ',trace[0][1])
    print()

    aper_flux_8 = psfcat_starlist['FLUX_APER'][:,7]
    index = np.argsort(aper_flux_8)
    # reverse order
    index = index[::-1]
    print('APER_FLUX_8: ', aper_flux_8[index[0]], aper_flux_8[index[-1]]) 
    print()

    psfcat_sample = psfcat_starlist[index[0:99]]

    print('data.meta: ', psfcat_sample.meta['INFILE'])

    key=raw_input("Enter any key to continue: ")

    #plot_vignette(data=psfcat_sample, lutscale='minmax') 
    plot_vignette(data=psfcat_sample) 
    print('plot vignette completed')

    key=raw_input("Enter any key to continue: ")
