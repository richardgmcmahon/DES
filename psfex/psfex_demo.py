"""

 Demo analysis DESDM psf products using Erin Sheldons psfex python function.

 See also:

 https://github.com/GalSim-developers/GalSim/blob/releases/1.1/galsim/des/des_psfex.py

 You will also need to be familiar with Emmanuel Bertin's psfex documentation
 http://www.astromatic.net/software/psfex

 rd_psf reads the DES products; currently the location
 and filename is hardwired

 Original version Richard McMahon (Nov 2013)

"""

# make future proof for Python 3
from __future__ import print_function, division

import os
import sys
import time

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
from astropy.table import Table
from astropy.io import fits
from astropy.stats import median_absolute_deviation

sys.path.append('/home/rgm/soft/python/lib/')
from psfex_util import *

#sys.path.append('/home/rgm/lib/python/psfex/')
import psfex

import galsim
print('galsim: ', galsim.__version__)
import galsim.des

mpl.use('Agg')

def psfex_galsim_compare(psfex_file=None, xpix=None, ypix=None):

  print('Reading: ', psfex_file)

  if xpix is None: xpix=5000.0
  if ypix is None: ypix=5000.0

  # make psf using galsim
  psfex_galsim=galsim.des.DES_PSFEx(psfex_file)
  xypos = galsim.PositionD(x=xpix, y=ypix)
  psf_galsim=psfex_galsim.getPSFArray(xypos)
  help(psf_galsim)

  #stamp_galsim = psf_galsim.drawImage()
  #help(stamp_galsim)

  # make psf using psfex
  pex = psfex.PSFEx(infile_psf)

  psf_psfex = pex.get_rec(xpix, ypix)

  plt.figure(figsize=(10.0, 8.0))
  plt.imshow(psf_psfex, interpolation='none')
  plt.colorbar()
  plt.ylabel('Pixels')
  plt.xlabel('Pixels')

  title= filename_image + ': ' + "{0:8.1f}".format(xpix) + "{0:8.1f}".format(ypix)
  title= RELEASE+ ': ' + filename_image 
  plt.title(title)
  #rgm.plotid()

  plotfile='psfex_demo_psfex_galsim_compare_psfex.png'
  print('Saving: ', plotfile)
  plt.savefig(plotfile)


  plt.figure(figsize=(10.0, 8.0))
  plt.imshow(psf_galsim, interpolation='none')
  plt.colorbar()
  plt.ylabel('Pixels')
  plt.xlabel('Pixels')

  title= filename_image + ': ' + "{0:8.1f}".format(xpix) + "{0:8.1f}".format(ypix)
  title= RELEASE+ ': ' + filename_image 
  plt.title(title)
  #rgm.plotid()

  plotfile='psfex_demo_psfex_galsim_compare_galsim.png'
  print('Saving: ', plotfile)
  plt.savefig(plotfile)


  plt.figure(figsize=(10.0, 8.0))  
  plt.imshow(psf_galsim-psf_psfex, interpolation='none')
  plt.colorbar()
  plt.ylabel('Pixels')
  plt.xlabel('Pixels')

  title= filename_image + ': ' + "{0:8.1f}".format(xpix) + "{0:8.1f}".format(ypix)
  title= RELEASE+ ': ' + filename_image 
  plt.title(title)
  #rgm.plotid()

  plotfile='psfex_demo_psfex_galsim_compare_difference.png'
  print('Saving: ', plotfile)
  plt.savefig(plotfile)


def rd_psf(infile=None, verbose=False, pause=False):
  """
  Read DES PSFEx psf file and do some diagnostic tests

  """

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

  # cycle through and plot each base function
  nbases=6
  for ibase in range(0,6): 
    #psf_mask.shape:  (1, 6, 25, 25)
    # slice out a single PCA component
    # The polynomial engine of PSFEx is the same as the one implemented 
    # the SCAMP software (Bertin 2006)
    image=psf_mask[0,ibase,:,:]
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

    plt.figure(figsize=(10.0, 8.0))

    lower=-1
    upper=10
    cmap='bone'
    plt.imshow(image, interpolation='none', cmap=cmap,
     vmin=data_min, vmax=data_max)

    #vmin=median+(lower*sigma_mad), vmax=median+(upper*sigma_mad))

    plt.colorbar()
    plt.ylabel('Pixels')
    plt.xlabel('Pixels')

    title= filename_image 

    title= RELEASE+ ': ' + filename_image + ' PSFEx base function: ' + str(ibase)

    plt.title(title)
    #rgm.plotid()


    plotfile='psfex_demo_base_function_'+ str(ibase) + '.png'
    print('Saving: ', plotfile)
    plt.savefig(plotfile)

    pex = psfex.PSFEx(infile_psf)

    return pex

  if pause: raw_input('Type any key to continue> ')

  # read the PSFEx .psf file
  pex = psfex.PSFEx(infile_psf)

  xpix=5000.0
  ypix=5000.0

  psf_pex = pex.get_rec(xpix, ypix)

  #print(psf_pex)


  # read the polynomial parameters
  #help(hdulist[1].header)
  polzero1=hdulist[1].header['POLZERO1']
  print('POLZERO1: ', polzero1)
  polzero2=hdulist[1].header['POLZERO2']
  print('POLZERO2: ', polzero2)
  polscal1=hdulist[1].header['POLSCAL1']
  print('POLSCAL1: ', polscal1)
  polscal2=hdulist[1].header['POLSCAL2']
  print('POLSCAL2: ', polscal2)

  xpsf=(xpix - polzero1)/polscal1
  ypsf=(ypix - polzero2)/polscal2
  print('xpix, ypix: ', xpix, ypix)
  print('xpsf, ypsf: ', xpsf, ypsf)

  # create model psf from base functions 
  # cst + x + x^2 + y + xy + y^2 
  psf = \
   psf_mask[0,0,:,:] + \
   (xpsf * psf_mask[0,1,:,:]) + \
   (xpsf*xpsf * psf_mask[0,2,:,:]) + \
   (ypsf * psf_mask[0,3,:,:]) + \
   (xpsf*ypsf * psf_mask[0,4,:,:]) + \
   (ypsf*ypsf * psf_mask[0,5,:,:])
     
  #print(psf)

  zoom_psf=True
  if zoom_psf:
    print('psf_sample: ', psf_sample)
    print('psf.shape: ', psf.shape)
    zoomed_psf = scipy.ndimage.interpolation.zoom(psf, psf_sample)
    print('zoomed_psf.shape: ', zoomed_psf.shape)
    #psf=zoomed_psf[1:-1,1:-1]
    psf=zoomed_psf
 
  plt.figure(figsize=(10.0, 8.0))

  plt.imshow(psf, interpolation='none')
  plt.colorbar()
  plt.ylabel('Pixels')
  plt.xlabel('Pixels')

  title= filename_image + ': ' + "{0:8.1f}".format(xpix) + "{0:8.1f}".format(ypix)
  title= RELEASE+ ': ' + filename_image 
  plt.title(title)
  #rgm.plotid()

  print('total flux in psf: ', np.sum(psf.flat))
  plotfile='psfex_demo_psf_1.png'
  print('Saving: ', plotfile)
  plt.savefig(plotfile)

  plt.figure(figsize=(10.0, 8.0))

  plt.imshow(psf_pex, interpolation='none')
  plt.colorbar()
  plt.ylabel('Pixels')
  plt.xlabel('Pixels')

  title= filename_image + ': ' + "{0:8.1f}".format(xpix) + "{0:8.1f}".format(ypix)
  title= RELEASE+ ': ' + filename_image 
  plt.title(title)
  #rgm.plotid()

  print('total flux in psf_pex: ', np.sum(psf_pex.flat))
  plotfile='psfex_demo_psf_2.png'
  print('Saving: ', plotfile)
  plt.savefig(plotfile)


  plt.figure(figsize=(10.0, 8.0))

  plt.imshow(psf-psf_pex, interpolation='none')
  plt.colorbar()
  plt.ylabel('Pixels')
  plt.xlabel('Pixels')

  title= filename_image + ': ' + "{0:8.1f}".format(xpix) + "{0:8.1f}".format(ypix)
  title= RELEASE+ ': ' + filename_image 
  plt.title(title)
  #rgm.plotid()

  plotfile='psfex_demo_psf_1_2.png'
  if zoom_psf: plotfile='psfex_demo_psf_1_2_zoomed.png'
  print('Saving: ', plotfile)
  plt.savefig(plotfile)

  data=psf

  plotfile='psfex_demo_radial_profile_psf.png'
  center=[12.0,12.0]
  radial_profile(data, center, sigma=None, plotfile=plotfile)

  data=psf_pex
  plotfile='psfex_demo_radial_profile_psf_pex.png'
  center=[12.0,12.0]
  radial_profile(data, center, sigma=None, plotfile=plotfile)

  return pex

def radial_profile(data, center, sigma=None, plotfile=None):
    """
    based on http://stackoverflow.com/questions/21242011/most-efficient-way-to-calculate-radial-profile

    """

    from scipy.stats import norm

    # there might be a numpy/scipy replacement of this
    import matplotlib.mlab as mlab

    y,x = np.indices((data.shape)) # first determine radii of all pixels
    print(min(x.flat), max(x.flat))
    print(min(y.flat), max(y.flat))
    r = np.sqrt((x-center[0])**2+(y-center[1])**2)    
    print(min(r.flat), max(r.flat))

    # radius of the image.
    r_max = np.max(r)  

    #ring_brightness, radius = np.histogram(r, weights=data, bins=r_max)
    #plt.plot(radius[1:], ring_brightness)

    # could sort on r and plot line
    print(len(r), r.shape)
    r=r.flat

    print(len(data), data.shape)
    data=data.flat

    index=np.argsort(r)
    print(len(index))
    r=r[index]
    print(len(r))
    data=data[index]
    print(len(data))

    plt.figure(figsize=(8.0, 8.0))

    plt.plot(r, data, 'k.')
    plt.xlim(0.0,12.5)
    plt.xlabel('radius (pixels)')  

    #gaussian = norm(loc = -1., scale = 1.0)
  
    pdf=False
    if pdf:
      xmin=0
      xmax=r_max

      nsteps=100
      x = np.linspace(xmin,xmax,nsteps+1)

      pdf=mlab.normpdf(x,0.0,sigma)
      xrange=xmax-xmin
      dx=xrange/nsteps
      cdf = np.cumsum(pdf*dx)
      plt.plot(x,cdf)

    if plotfile is None: plotfile='psfex_demo_radial_profile.png'
    print('Saving: ', plotfile)
    plt.savefig(plotfile)

    #plt.show()


def check_data(pause=True):

  # read the image and print some info
  infile=infile_image
  print('Read the image file: ', infile)
  hdulist = fits.open(infile)
  hdulist.info()

  if pause: raw_input('Type any key to continue> ')

  print(hdulist[0].header)
  print(hdulist[1].header)
  print(hdulist[2].header)

  if pause: raw_input('Type any key to continue> ')

  # read the psfcat file
  infile=infile_psfcat
  print('Read the psfcat file: ', infile)
  hdulist = fits.open(infile)
  hdulist.info()

  if pause: raw_input('Type any key to continue> ')

  print(hdulist[0].header)
  print(hdulist[1].header)
  print(hdulist[2].header)

  if pause: raw_input('Type any key to continue> ')

  # read the data 
  data=hdulist[2].data
  cols = data.columns
  print(cols.info)
  print(cols.names)
  infile=psf
  hdulist = fits.open(infile)
  hdulist.info()
  print(hdulist[0].header)
  print(hdulist[1].header)


  xdata=data['FLUX_APER']
  xdata=data['FLUX_RADIUS']

  # aper 3
  ydata=data['MAG_APER'][0:,2]

  #add a nominal zeropoint
  #ydata=ydata[0:,2]+25.0
  ydata=ydata+25.0

  print('len(xdata), len(ydata): ', len(xdata), len(ydata))
  print('min(xdata), max(xdata): ', min(xdata), max(xdata))
  print('min(ydata), max(ydata): ', min(ydata), max(ydata))


  #ydata_snr=data['FLUX_APER']
  #ydata_snr=ydata_snr[0:,2]

  plt.figure(figsize=(8.0, 8.0))

  #plt.plot(xdata, ydata, '.k')
  ndata=len(xdata)
  plt.scatter(xdata, ydata, marker='.', s=1, label=str(ndata))
  plt.xlim(-1.0,10.0)
  plt.ylim(20.0,5.0)

  plt.ylabel('MAG_APER_3 (uncalibrated)')
  plt.xlabel('FLUX_RADIUS (pixels)')
  plt.legend()


  title= RELEASE+ ': ' + filename_image 
  plt.title(title)

  #rgm.plotid()

  plotfile='psfex_demo_flux_radius_v_mag_aper.png'

  plt.savefig(plotfile)

def run_demo():

  print('run_demo')

  xcoords=[1000.0, 2000.0, 3000.0, 5000.0, 7000.0, 9000.0]
  ycoords=[1000.0, 2000.0, 3000.0, 5000.0, 7000.0, 9000.0]

  nx=9
  ny=9
  xcoords=np.linspace(1000.0, 9000.0, num=nx)
  ycoords=np.linspace(1000.0, 9000.0, num=ny)

  image=np.zeros((25*nx,25*ny), dtype=np.float64)
  print('image.shape: ', image.shape)
  print('image.type:  ', image.size)
  print('image.dtype: ', image.dtype)

  pex = psfex.PSFEx(psf)
  print('pex.get_center(0,0): ', pex.get_center(0,0))

  x0=-25
  for x in xcoords:
    y0=-25
    x0=x0+25
    for y in ycoords:
      y0=y0+25
      print(x, y, x0, y0)
      image[x0:x0+25,y0:y0+25] = pex.get_rec(x, y)
     
  plt.figure(figsize=(10.0, 8.0))

  print(sum(image.flat))

  plt.imshow(image, interpolation='none')
  plt.colorbar()
  plt.ylabel('Pixels')
  plt.xlabel('Pixels')

  title= filename_image + ': ' + "{0:8.1f}".format(x) + "{0:8.1f}".format(y)

  title= RELEASE+ ': ' + filename_image 

  plt.title(title)
  #rgm.plotid()


  plotfile='psfex_demo_image_psf.png'
  print('Saving: ', plotfile)
  plt.savefig(plotfile)

  image0= pex.get_rec(5000.0, 5000.0)
  image0 = np.tile(image0, (nx,ny))
  print(image0.shape)
  print(image0.size)
  print(image0.dtype)

  print('plot the psf residuals')
  imagediff=image-image0

  plt.figure(figsize=(10.0, 8.0))

  plt.imshow(imagediff, interpolation='none')
  plt.colorbar()
  plt.ylabel('Pixels')
  plt.xlabel('Pixels')

  title= filename_image + ': ' + "{0:8.1f}".format(x) + "{0:8.1f}".format(y)
  title= RELEASE+ ': ' + filename_image 
  plt.title(title)
  #rgm.plotid()

  plotfile='psfex_demo_image_psf_residuals.png'
  print('Saving: ', plotfile)
  plt.savefig(plotfile)

if __name__ == '__main__':
    """

    """

    t0=time.time()

    RELEASE='SVA1'

    TILE='DES1000+0209'
    datapath='/data/desardata/SVA1/COSMOS/' + TILE + '/'

    TILE='DES0449-4706'
    datapath='/data/desardata/SVA1/' + TILE + '/'
  

    filename_psf= TILE + '_i_psfcat.psf' 
    infile_psf=datapath + filename_psf

    filename_image= TILE+ '_i.fits.fz'
    infile_image= datapath + filename_image

    infile_psfcat= datapath + TILE + '_i_psfcat.fits.fz'
    psf= datapath + TILE + '_i_psfcat.psf'

    psfex_galsim_compare(psfex_file=infile_psf)

    rd_psf(infile_psf)

    print('Elapsed time(secs): ',time.time() - t0)

    pause=True
    if pause: raw_input('Type any key to continue> ')

    plot_vignette()
    print('Elapsed time(secs): ',time.time() - t0)

    check_data()
    print('Elapsed time(secs): ',time.time() - t0)

    run_demo()
    print('Elapsed time(secs): ',time.time() - t0)

    plot_radial_profile=True
    if plot_radial_profile:
        pex=rd_psf(infile_psf)
        print('fwhm: ', pex.get_fwhm())
        sigma=pex.get_fwhm()/2.3548
        psf=pex.get_rec(5000.0, 5000.0)
        print('psf.shape: ', psf.shape)
        #help(psf)

        center=[12.0,12.0]
        radial_profile(psf,center=center,sigma=sigma)

        debug=False
  

