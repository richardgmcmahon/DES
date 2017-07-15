"""
Simulate quasar pairs to investigate DESDM completeness

2 issue that we want to invesigate:


(i)    what is the classifier produce for close pairs
(ii)   what is the reliabiity of the photometry produced


Flow:

make to offset 25x25 psfs and place into a 60 x 60 array

for DES this is into a 15" x 15" array assuming 0.27" per pixel

TODO:

tidy up the broadcasting and other array indices

BEWARE off-by-one and counting issues

e.g. if we have the default psf 25x25 array and drop a pair
into a 50x50 array centred on one psf the offset can have a 
maximum of of 12 pixels 

range 0:24 with 12 as the centre (i.e 13th pixel)

adding an offset source with a 25x25 that is 13 pixels off will 
cover start on 25:49 hitting the edge; 14 pixels off will hit
the edge and generate a broadcasting error.


"""

import sys

import os
import sys

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

#sys.path.append('/home/rgm/soft/OM10/OM10/om10/')
#import om10

sys.path.append('/home/rgm/lib/python/psfex/')

import psfex
from astropy.io import fits

import matplotlib.pyplot as plt

import numpy as np

sys.path.append('/home/rgm/soft/python/lib/rgm/')
from plotid import *

def make_pair(pex=None, xcoord=0, ycoord=0, 
 xoffset=0, yoffset=0, xsize=25, ysize=25, ratio=1.0):
  """

   Assumes PSFEx default psf cutout size of 25x25. 

   This could be more flexible to take account of a range
   of offsets. e.g. the output size could be 25+xoffset, 25+yoffset
   in size rather than current default of 61x61 which allows
   a maximum of [x,y]offset of 25. 

   simulate a pair into a 61 x 61 pixel array
   make it odd so that centre is on a pixel centre.

   psfex makes 2d arrays 25x25 pixels

   initial version uses psfex at each point but it might be faster 
   to shift array by integer number of pixels

   also for generic pixel location use:
   scipy.ndimage.interpolation.shift with spline 3. 

   Lanczos interpolation would be ideal but the above should be 
   good enough. SWARP or CASU software could also be possible options.


  """

  PSFEx_PSF_SIZE=(25,25)
  # make it odd so that centre is on a pixel centre
  xsize_sim=61
  ysize_sim=61

  image1=np.zeros((xsize_sim,ysize_sim), dtype=np.float64)
  image2=np.zeros((xsize_sim,ysize_sim), dtype=np.float64)

  print xcoord, ycoord
  sim=pex.get_rec(xcoord, ycoord)
  print sim.shape
  print 'xoffset: ', xoffset  

  # assumes xcen, ycen = 30, 30 starting from 0,0 for 61x61
  # i.e. (size-1)/2
  # 25x25 inset is 12,12 -> 30,30
  #                 0,0  -> 18,18
  #                24,24 -> 42,42
  # note slicing is 18:43
  print   image1[18:43,18:43].shape
  image1[18:43,18:43] = sim
  print image2[18+xoffset:43+xoffset,18:43].shape
  image2[18+xoffset:43+xoffset,18:43] = sim

  sim=image1+image2

  return sim

def mymad(data, median=None, sigma=False):
  """
  compute median absolute deviation
  Options: 
    provide precomputed median
    return the equivalenet sigma
    maybe offer variance too
  """
  if median is None: median=np.median(data)

  mad=np.median (abs(data-median))

  if sigma: mad=mad/0.6745

  return mad


RELEASE='SVA1'
datapath='/data/desardata/SVA1/COSMOS/DES1000+0209/'

filename_psf='DES1000+0209_i_psfcat.psf' 
infile_psf=datapath + filename_psf

# read the psf file
infile=infile_psf
hdulist = fits.open(infile)
hdulist.info()
print 'len(hdulist): ', len(hdulist)
print hdulist[0].header
print hdulist[1].header


# read the image file
filename_image='DES1000+0209_i.fits.fz'
infile_image= datapath + filename_image
hdulist_image = fits.open(infile_image)
hdulist_image.info()
print 'len(hdulist_image): ', len(hdulist_image)
print hdulist_image[0].header
print hdulist_image[1].header
print hdulist_image[2].header

data=hdulist_image[1].data
print 'Computing Median  '

median=np.median(data.flat)
print
print 'Median:  ', median
print 'sigma(MAD):     ', mymad(data, median=median, sigma=True)
print 'sigma(MAD):     ', mymad(data.flat, median=median, sigma=True)
print 'Minimum: ', np.min(data)
print 'Maximum: ', np.max(data)


data=hdulist_image[2].data
median=np.median(data.flat)
print
print 'Median: ', median
print 'sigma(MAD):     ', mymad(data, median=median, sigma=True)
print 'sigma(MAD):     ', mymad(data.flat, median=median, sigma=True)
print 'Minimum: ', np.min(data)
print 'Maximum: ', np.max(data)

title= RELEASE+ ': ' + filename_psf

xcoord=5000.0
ycoord=5000.0

# DES COADD image pixel size is 0.27 arc seconds per pixel

offsets=np.linspace(0.1, 2.0, num=20)

# make image
nx=5
ny=5
# 4 times 61 plus a margin
xsize=250
ysize=250

#xsize=1000
#ysize=1000

bigimage=np.zeros((xsize,ysize), dtype=np.float64)
print bigimage.shape
print bigimage.size
print bigimage.dtype

pex = psfex.PSFEx(infile_psf)

# range from 0 to 16 pixels i.e 0 to 4.32 assuming 0.27 arc sec per pixel 
xoffsets=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]

xoffset=5
# make a pair as a test
sim=make_pair(pex=pex, xcoord=xcoord, ycoord=ycoord, 
 xoffset=xoffset, xsize=xsize, ysize=ysize, ratio=1.0)
print sim.shape
print sim.size
print 'total flux in sim: ', np.sum(sim)

x0=-25


#
xstep=61
ystep=61

#xstep=200
#ystep=200

# create mosaic of pair sims
ioffset=-1
for ix in range(0, 4):
  for iy in range(0, 4):
    ioffset=ioffset+1
    xoffset=xoffsets[ioffset]
    print 'ioffset, ix, iy: ', ioffset, ix, iy, xoffset

    sim=make_pair(pex=pex, xcoord=xcoord, ycoord=ycoord, 
     xoffset=xoffset, xsize=xsize, ysize=ysize, ratio=1.0)
    
    ixstart=ix*xstep
    ixend=ixstart+xstep
    iystart=iy*ystep
    iyend=iystart+ystep
  
    print 'ixstart, ixend: ', ixstart, ixend
    print 'iystart, iyend: ', iystart, iyend
    print sim.shape
    print sim.size

    print 'bigimage.shape: ', bigimage.shape
    print 'bigimage.size: ', bigimage.size

    print bigimage[ixstart:ixend:1, iystart:iyend:1].shape
    print bigimage[ixstart:ixend:1, iystart:iyend:1].size

    # paste sim into image
    bigimage[ixstart:ixend:1, iystart:iyend:1] = sim

# save the 250x250 simulation with 16 pairs
image=bigimage

# replicate 250x250 into a 10,000 by 10,000 image i.e. 40 x 40
xsize=10000
ysize=10000
sim_image=np.zeros((xsize,ysize), dtype=np.float64)
# create mosaic of sims
ioffset=-1
xstep=250
ystep=250
nsteps=4
for ix in range(0, nsteps):
  for iy in range(0, nsteps):

    ixstart=ix*xstep
    ixend=ixstart+xstep

    iystart=iy*ystep
    iyend=iystart+ystep
  
    print 'ixstart, ixend: ', ixstart, ixend
    print 'iystart, iyend: ', iystart, iyend

    sim_image[ixstart:ixend:1, iystart:iyend:1] = image

image=sim_image
  
plt.figure(figsize=(10.0, 8.0))

print sum(image.flat)

plt.imshow(image, interpolation='none')
plt.colorbar()
plt.ylabel('Pixels')
plt.xlabel('Pixels')
plotid(progname=True)

title= RELEASE+ ': ' + filename_psf

plt.title(title)
#rgm.plotid()

plotfile='psfex_sims_1.png'

plt.savefig(plotfile)


# scale the sim to have a mag of 20 using zeropoint;
# sims are normalised to have a total flux of 1.0
  
# write the model to a fits file with the DES header


# read the image file
filename_out='image_sim.fits'
hdu=fits.PrimaryHDU(image)
hdulist=fits.HDUList([hdu])
hdulist.writeto(filename_out)

