"""

  some utility functions for psfs based on PSFEx. Initially focus is on
  DES data but also VISTA, VST, Gaia, LSST will be considered

  see also psfex_demo

  could read a PSFex parameter file


"""
from __future__ import print_function, division

import os
import sys

import matplotlib.pyplot as plt

import numpy as np

import psfex

from skimage import measure

from astropy.table import Table, join
from astropy.io import fits
from astropy.stats import median_absolute_deviation

sys.path.append('/home/rgm/soft/python/lib/')
from librgm.plotid import plotid
#from plotid import *

def set_des_apertures(unit='pixels', verbose=False):
    """

    see https://opensource.ncsa.illinois.edu/confluence/display/DESDM/Overview+of+data+columns+in+DES+Databases#OverviewofdatacolumnsinDESDatabases-MAG_APER#_GRIZY

    0.270 arcsec/pixel
    1.00  arcsec = 3.704 pixels
 
          Diameter   Diameter
    Aperture  [arcsec]   [pixels]
   Mag_aper1    0.5        1.85
   Mag_aper2    1.0        3.70
   Mag_aper3    1.5        5.55
   Mag_aper4    2.0        7.41
   Mag_aper5    3.0       11.11
   Mag_aper6    4.0       14.81
   Mag_aper7    5.0       18.52
   Mag_aper8    6.0       22.22
   Mag_aper9    7.0       25.93
   Mag_aper10   8.0       29.63
   mag_aper11  12.0       44.44
   mag_aper12  18.0       66.67

    """

    aperture_diameter_pixels = np.asarray([
        1.85,
        3.70, 
        5.55, 
        7.41, 
        11.11, 
        14.81, 
        18.52,
        22.22,
        25.93,
        29.63,
        44.44,
        66.67 ])

    aperture_diameter_arcsecs = np.asarray([
         0.50,
         1.00,
         1.50,
         2.00,
         3.00,
         4.00,
         5.00,
         6.00,
         7.00,
         8.00,
        12.00,
        18.00])

    if verbose:
        np.set_printoptions(precision=3)
        print(aperture_diameter_arcsecs/aperture_diameter_pixels)

        np.set_printoptions(precision=3)
        print(aperture_diameter_pixels/aperture_diameter_arcsecs)

    aperture_diameter = aperture_diameter_pixels

    if unit != 'pixels':     aperture_diameter = aperture_diameter_arcsecs

    return  aperture_diameter

def rd_psf(infile=None, plots=False,
    prefix=None, suffix=None, verbose=False, pause=False):
    """
    Read DES PSFEx psf file and do some diagnostic tests

    """

    if prefix is None: prefix = ""

    # read the psf file
    hdulist = fits.open(infile)
    hdulist.info()
    print(hdulist[0].header)
    print(hdulist[1].header)

    data=hdulist[1].data
    psf_mask=data['PSF_MASK']

    print('len(psf_mask.flat): ', len(psf_mask.flat))
    print('psf_mask.shape: ', psf_mask.shape)
    print('psf_mask.type:  ', psf_mask.size)
    print('psf_mask.dtype: ', psf_mask.dtype)

    print('min: ', min(psf_mask.flat))
    print('max: ', max(psf_mask.flat))
    print('len: ', len(psf_mask.flat))

    psf_sample = hdulist[1].header["PSF_SAMP"]

    # cycle through and plot each base functions
    
    if plots: 
        plot_psfex_basefunctions(psf_bases=psf_mask, infile=infile,
            prefix=prefix, pause=pause)

    pex = psfex.PSFEx(infile)

    print('fwhm: ', pex.get_fwhm())
    sigma=pex.get_fwhm()/2.3548
    psf=pex.get_rec(0.0, 0.0)
    print('psf.shape: ', psf.shape)

    return pex


def rd_psfex_starlist(infile=None, pause=False):
    """

    """


    # read the psfcat file
    print('Read the psfex starlist file: ', infile)
    hdulist = fits.open(infile)
    print('Number of extensions: ', len(hdulist))
    hdulist.info()

    print(hdulist[0].header)
    print(hdulist[1].header)
    print(hdulist[2].header)

    # read the data 
    data=hdulist[2].data
    cols = data.columns
    print(cols.info)
    print(cols.names)

    print('Number of rows: ', len(data))

    table = Table.read(infile, hdu=2)    

    table.meta['INFILE']=infile

    print( 'table metadata: ', table.meta)

    if pause: raw_input('Type any key to continue> ')

    return table

def rd_psfcat(infile=None, pause=False, plots=True):
    """


    """
    # read the psfcat file
    print('Read the psfcat file: ', infile)
    hdulist = fits.open(infile)
    print('Number of extensions: ', len(hdulist))
    hdulist.info()

    print(hdulist[0].header)
    print(hdulist[1].header)
    print(hdulist[2].header)

    # read the data 
    data=hdulist[2].data
    cols = data.columns
    print(cols.info)
    print(cols.names)

    table = Table.read(infile, hdu=2)    

    table.meta['INFILE']=infile

    print( 'table metadata: ', table.meta)

    print('Number of rows: ', len(data))

    if pause: raw_input('Type any key to continue> ')

    if plots: 
        psfcat_plots(data=table, infile=infile)
        if pause: raw_input('Type any key to continue> ')

    return table


def psfcat_plots(data=None, infile=None, label=None, overplot=False,
    FLUX_SAT=None, SAMPLE_MINSN=None, SAMPLE_MAXELLIP=None,
    RANGE_FLUX_RADIUS=None,
    RANGE_FLUX_APER=None,
    RANGE_ELONGATION=None):
    """

    could try to get the filename from the table metadata


    """

    try:
        print('data.meta: ', data.meta['INFILE'])
        print('overplot: ', overplot)
        raw_input('Type any key to continue> ')
    except:
        pass

    filename = os.path.basename(infile)  

    plotfile_suffix = ''
    if label is not None: plotfile_suffix = '_'+ label
    
    # PSFEx sample selection config parameters
    if FLUX_SAT is None: FLUX_SAT = 65000
    if SAMPLE_MINSN is None: SAMPLE_MINSN = 20
    if SAMPLE_MAXELLIP is None: SAMPLE_MAXELLIP = 0.3

    # give the plot a unique figname thing
    #if not overplot: plt.clf("fig_psfcat_1")
    plt.figure("fig_psfcat_1", figsize=(8.0, 8.0))

    xdata=data['FLUX_RADIUS']

    # extract a flux/mag vector 
    ydata=data['FLUX_APER'][:,7]

    #ydata=data['MAG_APER'][:,3]

    print('min(ydata), max(ydata): ', min(ydata), max(ydata))

    #add a nominal zeropoint
    #ydata=ydata[0:,2]+25.0
    ydata = np.log10(ydata)

    print('len(xdata), len(ydata): ', len(xdata), len(ydata))
    print('min(xdata), max(xdata): ', min(xdata), max(xdata))
    print('min(ydata), max(ydata): ', min(ydata), max(ydata))


    #plt.plot(xdata, ydata, '.k')
    ndata=len(xdata)

    if not overplot: 
        plt.scatter(xdata, ydata, marker='.', s=1, color='red',
            label=str(ndata))

    if overplot: 
        plt.scatter(xdata, ydata, marker='.', s=1, color='green',
            label=str(ndata))

    xrange=(-1.0,30.0)
    if RANGE_FLUX_RADIUS is not None: xrange = RANGE_FLUX_RADIUS

    yrange=(-1.0,8.0)
    if RANGE_FLUX_APER is not None: yrange = np.log10(RANGE_FLUX_APER)

    plt.xlim(xrange)
    plt.ylim(yrange)

    #xline=[np.log10(2.0), np.log10(2.0)]
    #yline=yrange
    #plt.plot(xline, yline, color='red', linestyle='--')

    plt.ylabel('Log(FLUX_APER_8) (uncalibrated)')
    plt.xlabel('FLUX_RADIUS_HALF_LIGHT (pixels)')
    #plt.grid(which='both')
    plt.grid(True)

    plt.legend()
    title= infile
    plt.title(title, fontsize='medium')
    plotid(progname=True)

    plotfile = filename + '_flux_radius_v_flux_aper_8' + \
        plotfile_suffix + '.png'

    plt.savefig(plotfile)
    #plt.clf()
    print('Saving: ', plotfile)

   
    # need to look also at: SNR_WIN
    # could also look at aper 4 
    # give the plot a unique figname thing
    #if not overplot: plt.clf("fig_psfcat_2")
    plt.figure("fig_psfcat_2", figsize=(8.0, 8.0))

    # extract a flux/mag vector for aper 8
    xdata=data['FLUX_APER'][:,7]/data['FLUXERR_APER'][:,7]
    ydata=data['FLUX_APER'][:,7]

    print('min(ydata), max(ydata): ', min(ydata),max(ydata))

    #add a nominal zeropoint
    #ydata=ydata[0:,2]+25.0
    xdata = np.log10(xdata)
    ydata = np.log10(ydata)

    print('len(xdata), len(ydata): ', len(xdata), len(ydata))
    print('min(xdata), max(xdata): ', min(xdata), max(xdata))
    print('min(ydata), max(ydata): ', min(ydata), max(ydata))

    #plt.plot(xdata, ydata, '.k')
    ndata=len(xdata)

    ndata=len(xdata)

    if not overplot: 
        plt.scatter(xdata, ydata, marker='.', s=1, color='red',
            label=str(ndata))

    if overplot: 
        plt.scatter(xdata, ydata, marker='.', s=1, color='green',
            label=str(ndata))

    xrange=(-2.0,6.0)
    yrange=(-1.0,8.0)

    plt.xlim(xrange)

    if RANGE_FLUX_APER is not None: yrange = np.log10(RANGE_FLUX_APER)
    plt.ylim(yrange)

    xline = [np.log10(SAMPLE_MINSN), np.log10(SAMPLE_MINSN)]
    yline = yrange
    plt.plot(xline, yline, color='red', linestyle='--')

    plt.xlabel('S/N')
    plt.ylabel('Log(FLUX_APER_8) (uncalibrated)')
    plt.legend()
    title= infile
    plt.title(title, fontsize='medium')
    plotid(progname=True)

    plotfile= filename + '_fluxSN_v_flux_aper_8' + \
        plotfile_suffix + '.png'

    plt.savefig(plotfile)
    #plt.clf()
    print('Saving: ', plotfile)


    # 
    #if not overplot: plt.clf()
    plt.figure("fig_psfcat_3", figsize=(8.0, 8.0))

    print('Plot FLUX_MAX versus FLUX')

    # extract a flux  vector 
    xdata=data['FLUX_MAX']
    ydata=data['FLUX_APER'][:,7]

    print('len(xdata), len(ydata): ', len(xdata), len(ydata))
    print('min(xdata), max(xdata): ', min(xdata), max(xdata))
    print('min(ydata), max(ydata): ', min(ydata), max(ydata))

    xdata=np.log10(xdata)
    ydata=np.log10(ydata)

    print('len(xdata), len(ydata): ', len(xdata), len(ydata))
    print('min(xdata), max(xdata): ', min(xdata), max(xdata))
    print('min(ydata), max(ydata): ', min(ydata), max(ydata))

    if not overplot: 
        plt.scatter(xdata, ydata, marker='.', s=1, color='red',
            label=str(ndata))

    if overplot: 
        plt.scatter(xdata, ydata, marker='.', s=1, color='green',
            label=str(ndata))

    plt.xlabel('Log(FLUX_MAX) (uncalibrated)')
    plt.ylabel('Log(FLUX_APER_8) (uncalibrated)')
    plt.legend()

    yline=yrange
    xline=[np.log10(FLUX_SAT), np.log10(FLUX_SAT)]
    plt.plot(xline, yline, color='red', linestyle='--')

    xline=xrange
    yline=[np.log10(FLUX_SAT), np.log10(FLUX_SAT)]
    plt.plot(xline, yline, color='red', linestyle='--')


    title= infile
    plt.title(title, fontsize='medium')
    plotid(progname=True)

    plotfile= filename + '_flux_max_v_flux_aper_8' + \
        plotfile_suffix + '.png'

    plt.savefig(plotfile)
    #plt.clf()
    print('Saving: ', plotfile)

    #if not overplot: plt.clf()
    plt.figure("fig_psfcat_4", figsize=(8.0, 8.0))

    print('Plot ELONGATION versus FLUX')

    # extract a flux  vector 
    xdata=data['ELONGATION']
    ydata=data['FLUX_APER'][:,7]

    print('len(xdata), len(ydata): ', len(xdata), len(ydata))
    print('min(xdata), max(xdata): ', min(xdata), max(xdata))
    print('min(ydata), max(ydata): ', min(ydata), max(ydata))
    print()
    ydata=np.log10(ydata)

    print('len(xdata), len(ydata): ', len(xdata), len(ydata))
    print('min(xdata), max(xdata): ', min(xdata), max(xdata))
    print('min(ydata), max(ydata): ', min(ydata), max(ydata))
    print()

    if not overplot: 
        plt.scatter(xdata, ydata, marker='.', s=1, color='red',
            label=str(ndata))

    if overplot: 
        plt.scatter(xdata, ydata, marker='.', s=1, color='green',
            label=str(ndata))

    plt.xlabel('ELONGATION')
    plt.ylabel('Log(FLUX_APER_8) (uncalibrated)')

    xrange=(0.0,10.0)
    yrange=(-1.0,8.0)

    if RANGE_FLUX_APER is not None: yrange = np.log10(RANGE_FLUX_APER)
    if RANGE_ELONGATION is not None: xrange = RANGE_ELONGATION

    plt.xlim(xrange)
    plt.ylim(yrange)

    plt.grid(True)

    plt.legend()

    #xline=xrange
    #yline=[np.log10(FLUX_SAT), np.log10(FLUX_SAT)]
    #plt.plot(xline, yline, color='red')


    title= infile
    plt.title(title, fontsize='medium')
    plotid(progname=True)

    plotfile= filename + '_elongation_v_flux_aper_8' + \
        plotfile_suffix + '.png'

    plt.savefig(plotfile)
    #plt.clf()
    print('Saving: ', plotfile)


    # plot SNR_WIN versus FLUX

    #if not overplot: plt.clf()
    plt.figure("fig_psfcat_5", figsize=(8.0, 8.0))

    print('Plot SNR_WIN versus FLUX')

    # extract a flux  vector 
    xdata=data['SNR_WIN']
    ydata=data['FLUX_APER'][:,7]

    print('len(xdata), len(ydata): ', len(xdata), len(ydata))
    print('min(xdata), max(xdata): ', min(xdata), max(xdata))
    print('min(ydata), max(ydata): ', min(ydata), max(ydata))
    print()
    xdata=np.log10(xdata)
    ydata=np.log10(ydata)

    print('len(xdata), len(ydata): ', len(xdata), len(ydata))
    print('min(xdata), max(xdata): ', min(xdata), max(xdata))
    print('min(ydata), max(ydata): ', min(ydata), max(ydata))
    print()

    if not overplot: 
        plt.scatter(xdata, ydata, marker='.', s=1, color='red',
            label=str(ndata))

    if overplot: 
        plt.scatter(xdata, ydata, marker='.', s=1, color='green',
            label=str(ndata))

    plt.xlabel('Log(SNR_WIN)')
    plt.ylabel('Log(FLUX_APER_8) (uncalibrated)')

    xrange=(-2.0,6.0)
    yrange=(-1.0,8.0)
    if RANGE_FLUX_APER is not None: yrange = np.log10(RANGE_FLUX_APER)

    plt.xlim(xrange)
    plt.ylim(yrange)

    plt.grid(True)

    plt.legend()

    xline = [np.log10(SAMPLE_MINSN), np.log10(SAMPLE_MINSN)]
    yline = yrange
    plt.plot(xline, yline, color='red', linestyle='--')

    #xline=xrange
    #yline=[np.log10(FLUX_SAT), np.log10(FLUX_SAT)]
    #plt.plot(xline, yline, color='red')


    title= infile
    plt.title(title, fontsize='medium')
    plotid(progname=True)

    plotfile= filename + '_SNR_WIN_v_flux_aper_8' + \
        plotfile_suffix + '.png'

    plt.savefig(plotfile)
    #plt.clf()
    print('Saving: ', plotfile)



    # plot SNR_WIN versus FLUX

    #if not overplot: plt.clf()
    plt.figure("fig_psfcat_6", figsize=(8.0, 8.0))

    print('Plot FLUX_APER/FLUX_MAX versus FLUX')

    # extract a flux  vector 
    xdata=data['FLUX_APER'][:,7]/data['FLUX_MAX']
    ydata=data['FLUX_APER'][:,7]


    print('len(xdata), len(ydata): ', len(xdata), len(ydata))
    print('min(xdata), max(xdata): ', min(xdata), max(xdata))
    print('min(ydata), max(ydata): ', min(ydata), max(ydata))
    print()
    xdata=np.log10(xdata)
    ydata=np.log10(ydata)

    print('len(xdata), len(ydata): ', len(xdata), len(ydata))
    print('min(xdata), max(xdata): ', min(xdata), max(xdata))
    print('min(ydata), max(ydata): ', min(ydata), max(ydata))
    print()

    if not overplot: 
        plt.scatter(xdata, ydata, marker='.', s=1, color='red',
            label=str(ndata))

    if overplot: 
        plt.scatter(xdata, ydata, marker='.', s=1, color='green',
            label=str(ndata))

    plt.xlabel('Log(FLUX_APER_8/FLUX_MAX)')
    plt.ylabel('Log(FLUX_APER8) (uncalibrated)')

    xrange=( 0.5,2.5)
    yrange=(-1.0,8.0)
    if RANGE_FLUX_APER is not None: yrange = np.log10(RANGE_FLUX_APER)

    plt.xlim(xrange)
    plt.ylim(yrange)

    plt.grid(True)

    plt.legend()

    #xline=xrange
    #yline=[np.log10(FLUX_SAT), np.log10(FLUX_SAT)]
    #plt.plot(xline, yline, color='red')


    title= infile
    plt.title(title, fontsize='medium')
    plotid(progname=True)

    plotfile= filename + '_MAX_FLUX_over_TOTAL_FLUX_v_FLUX_APER_8' + \
        plotfile_suffix + '.png'

    plt.savefig(plotfile)
    #plt.clf()
    print('Saving: ', plotfile)




def plot_psfex_basefunctions(psf_bases=None, infile=None, 
    pause=False, prefix=None):
    """plot each psfex base functions

    """

    import numpy as np

    nbases=6
    for ibase in range(0,6): 
        #psf_mask.shape:  (1, 6, 25, 25)
        # slice out a single PCA component
        # The polynomial engine of PSFEx is the same as the one implemented 
        # the SCAMP software (Bertin 2006)
        image=psf_bases[0,ibase,:,:]
        if pause: raw_input('Type any key to continue> ')

        print('min, max: ', min(image.flat), max(image.flat))

        data_min=min(image.flat)
        data_max=max(image.flat)

        median=np.median(image.flat)
        print('median: ', median)

        MAD=median_absolute_deviation(image.flat)
        print('MAD: ', MAD)
        sigma_mad=MAD*1.4826
        print('sigma_mad: ', sigma_mad)

        print('len(image.flat): ', len(image.flat))
        print('image.shape: ', image.shape)
        print('image.type:  ', image.size)
        print('image.dtype: ', image.dtype)

        plt.figure("rd_psf",figsize=(10.0, 8.0))

        lower=-1
        upper=10
        cmap='bone'
        plt.imshow(image, interpolation='none', cmap=cmap,
         vmin=data_min, vmax=data_max)

        #vmin=median+(lower*sigma_mad), vmax=median+(upper*sigma_mad))

        plt.colorbar()
        plt.ylabel('Pixels')
        plt.xlabel('Pixels')


        filename = os.path.basename(infile)  
        title= filename

        title= prefix + ': ' + filename + \
            ' PSFEx base function: ' + str(ibase)

        plt.title(title)
        #rgm.plotid()

        suffix=filename
        if suffix is None: suffix = ''
        plotfile='psfex_base_function_'+ str(ibase) + '_' + suffix + '.png'
        print('Saving: ', plotfile)
        plt.savefig(plotfile)

        plt.clf()

def plot_vignette(data=None, infile_psfcat=None, lutscale=None,
    pause=True, verbose=False, RELEASE=None):

    """

    PSFEx does not work directly on images. Instead, it operates on SExtractor 
    catalogues that have a small image ("vignette") recorded for each detection.

    """

    print('Running plot_vignette')
    if data is not None: infile_psfcat = data.meta['INFILE']

    print('infile_psfcat: ', infile_psfcat)

    filename_psfcat = os.path.basename(infile_psfcat)  

    if data is None and infile_psfcat is not None:
    
        # read the psfcat file
        infile=infile_psfcat
        print('Read the psfcat file: ', infile)
        hdulist = fits.open(infile)
        hdulist.info()

        if pause: raw_input('Type any key to continue> ')

        if verbose:
            print(hdulist[0].header)
            print(hdulist[1].header)
            print(hdulist[2].header)
            if pause: raw_input('Type any key to continue> ')

        # read the data 
        data=hdulist[2].data
        cols = data.columns
        print(cols.info)
        print(cols.names)

    vignettes=data['VIGNET']

    print(vignettes[0:0])
    print(len(vignettes))

    flux_aper=data['FLUX_APER'][0:,2]
    isort=np.argsort(flux_aper)

    # choose brightest
    ivignette=isort[-1]

    print('Flux_Aper: ', data['FLUX_APER'][ivignette])
    print('Flux radius: ',data['FLUX_RADIUS'][ivignette])
    image=vignettes[ivignette]

    print('min, max: ', min(image.flat), max(image.flat))
    data_max=max(image.flat)

    itest=image.flat < -1e10
    print ('min, max: ', min(image.flat[itest]), max(image.flat[itest]), 
        len(image.flat[itest]))

    itest=image.flat > -1e10
    print ('min, max: ', min(image.flat[itest]), max(image.flat[itest]), 
        len(image.flat[itest]))


    itest=image.flat > -1e10
    median=np.median(image.flat[itest])
    print('median: ', median)

    MAD=median_absolute_deviation(image.flat[itest])
    print('MAD: ', MAD)
    sigma_mad=MAD*1.4826
    print('sigma_mad: ', sigma_mad)

    print('len(image.flat): ', len(image.flat))
    print('type(image): ', type(image))
    print('image.shape: ', image.shape)
    print('image.size:  ', image.size)
    print('image.dtype: ', image.dtype)

    # convert to avoid floating point overflows
    image = image.astype(np.float64, copy=False)
    print('image.dtype: ', image.dtype)
    image.clip(0.0)

    # location of maximum peak value
    imax = np.argmax(image)
    xymax = np.unravel_index(imax, image.shape)
    print('Location of maximum value: ',xymax)
    ImageMax= np.max(image)
    print('peak value: ', ImageMax)
    itest = image > -10000.0
    print('total intensity: ', np.sum(image[itest]))
    itest = image < -10000.0
    image[itest] = 0.0
    print('total intensity: ', np.sum(image))

    print('min, max: ', np.min(image.flat), np.max(image.flat))

    moments = measure.moments(image, order=2)
    print()
    print('measure.moments: \n', moments)
    print()
    total = moments[0, 0]
    cr = moments[0, 1] / moments[0, 0]
    cc = moments[1, 0] / moments[0, 0]
    print('Total: ', total)
    print('Weighted centroid (cr, cc): ', cr, cc)
    print()

    moments_central = measure.moments_central(image, cr, cc, order=2)
    print('measure.moments_central: \n', moments_central)
    print()

    print('measure.moments_central/total: \n', moments_central/total)
    print()

    moments_normalized = measure.moments_normalized(moments_central, order=2)
    print('measure.moments_normalized: \n', moments_normalized)

    plt.figure(figsize=(10.0, 8.0))

    lower=-1
    upper=100

    if lutscale is 'minmax':
        plt.imshow(image, interpolation='nearest',
            vmin=median+(lower*sigma_mad), vmax=ImageMax)

    if lutscale is None:
        plt.imshow(image, interpolation='nearest',
            vmin=median+(lower*sigma_mad), vmax=median+(upper*sigma_mad))

    plt.colorbar()
    plt.ylabel('Pixels')
    plt.xlabel('Pixels')

    title= filename_psfcat

    title = infile_psfcat
    if RELEASE is not None: title= RELEASE + ': ' + filename_psfcat 

    plt.title(title)
    plotid(progname=True)

    plotfile = filename_psfcat + '_psfex_demo_vignette.png'
    print('Saving: ', plotfile)
    plt.savefig(plotfile)


def match_starlist(psfcat=None, starlist=None):
    """
 
    """
    # now match with the stars and filter the psfcat data and redo the plots
    # SOURCE_NUMBER
    # subtract 1 since psfcat is indexed by row starting at 0
    index = starlist['SOURCE_NUMBER'] - 1
    X_IMAGE_starlist = starlist['X_IMAGE']
    Y_IMAGE_starlist = starlist['Y_IMAGE']

    X_IMAGE_psfcat = psfcat['X_IMAGE'][index]
    Y_IMAGE_psfcat = psfcat['Y_IMAGE'][index]
   
    delta_X_IMAGE = X_IMAGE_psfcat - X_IMAGE_starlist 
    delta_Y_IMAGE = X_IMAGE_psfcat - X_IMAGE_starlist 

    xdata = delta_X_IMAGE
    ydata = delta_Y_IMAGE
    print('len(xdata), len(ydata): ', len(xdata), len(ydata))
    print('min(xdata), max(xdata): ', min(xdata), max(xdata))
    print('min(ydata), max(ydata): ', min(ydata), max(ydata))



    return psfcat[index]
