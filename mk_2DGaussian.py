import sys
import numpy as np

import matplotlib.pyplot as plt

sys.path.append('/home/rgm/soft/python/lib/')
from librgm.plotid import plotid

#
#from astropy.modeling import models 
#g2d = models.Gaussian2D(amplitude=1.0, x_mean=0.0, y_mean=0.0, 
# x_stddev=0.5, y_stddev=0.5)
#print(g2d)

def makeGaussian(size=3, fwhm = 3, center=None):
    """ Make a square gaussian kernel.

    size is the length of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    """

    x = np.arange(0, size, 1, float)
    y = x[:,np.newaxis]

    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]

    return np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)


def mk_textfile(data=None, size=1.0, fwhm=1.0):
  """

  """

  # add Gaussian specs to filename
  outfilename="mygauss_" + str(fwhm) + "_"+ str(size)+"x"+str(size) + \
   "_conv_" + datestamp + ".txt"

  header="CONV NORM\n" + \
    '# '+str(size)+'x'+str(size)+' convolution mask of a gaussian PSF ' + \
    'with FWHM = ' + str(fwhm) + ' pixels'
  np.savetxt(outfilename, data, delimiter=' ', fmt='%8.5f',
    header=header, comments='')  

def plot_2DGaussian(data):
   """

   """   

   for irow in range(size):
      plt.plot(data[irow,:])


   plt.title('mk_2DGaussian; slices in Y')
   plt.xlabel('X pixel')
   plt.ylabel('Value')

   plotfile="mygauss_" + str(fwhm) + "_"+ str(size)+"x"+str(size) + \
     "_conv_" + datestamp + ".png"
   plotid(progname=True)

   print 'Saving: ', plotfile
   plt.savefig(plotfile)

   #plt.show()
   plt.close()


if __name__ == "__main__":

  import time

  now = time.localtime(time.time()) 
  datestamp = time.strftime("%Y%m%d",now)

  size=3
  fwhm=2
  data=makeGaussian(size=size, fwhm = fwhm, center=None)
  print(data)
  mk_textfile(data=data, size=size, fwhm=fwhm)

  size=9
  fwhm=5
  data=makeGaussian(size=size, fwhm = fwhm, center=None)
  print(data)
  mk_textfile(data=data, size=size, fwhm=fwhm)


  # CONV NORM
  # 15x15 convolution mask of a gaussian PSF with FWHM = 8.0 pixels.
  size=15
  fwhm=8
  data=makeGaussian(size=size, fwhm = fwhm, center=None)
  mk_textfile(data=data, size=size, fwhm=fwhm)
  print('shape: ',data.shape)
  print('data.size: ', data.size)
  print('np.sum(data.flatten): ', np.sum(data.flatten))
  mkplot=True
  if mkplot: plot_2DGaussian(data)



  # CONV NORM
  # 15x15 convolution mask of a gaussian PSF with FWHM = 8.0 pixels.
  size=25
  fwhm=12
  data=makeGaussian(size=size, fwhm = fwhm, center=None)
  mk_textfile(data=data, size=size, fwhm=fwhm)
  print('shape: ',data.shape)
  print('data.size: ', data.size)
  print('np.sum(data.flatten): ', np.sum(data.flatten))
  mkplot=True
  if mkplot: plot_2DGaussian(data)

  key=raw_input("Enter any key to continue: ")


  size=9
  fwhm=5
  data=makeGaussian(size, fwhm = fwhm, center=None)
  print(data)

  size=7
  fwhm=4
  data=makeGaussian(size, fwhm = fwhm, center=None)
  print(data)



