"""

WARNING: This currently does not work near 24hrs

RADEC_SOURCE_BAND in COADD table specifies the band used for RA, Dec in
table; need to check this.

also look at desdb for downloading on the fly option

see http://docs.astropy.org/en/stable/coordinates/matchsep.htm

TODO: go from COADD source to Single Epoch source and the
chip it FINALCUT Image and Chip it comes from.

https://opensource.ncsa.illinois.edu/confluence/display/DESDM/Y1A1+release+notes


Table Y1A1_COADD has a COADD_OBJECT
Table Y1A1_OBJECTS is the Main Single Epoch quantities table
(Y1A1_FINALCUT is a  synonym for Y1A1_OBJECTS)

DESDM data:

Y1A1: SE

images in OPS/red/'run'/red
segmentation images are in QA

RA, Dec analysis

could also compare the WIN and non-WIN positions

ALPHAWIN_J2000_[GRIZY], DELTAWIN_J2000_[GRIZY]

RADEC_SOURCE_BAND



"""
from __future__ import (print_function, division)

import os
import sys
import time

import numpy as np
from numpy import ma

import matplotlib as mpl
from matplotlib import pyplot as plt

from matplotlib import patches
from matplotlib.patches import Circle, Ellipse

from astropy.io import fits
from astropy.table import Table

from astropy import units as u
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS

from astropy import wcs

from astropy.stats import median_absolute_deviation
from astropy.stats import mad_std

sys.path.append('/home/rgm/soft/python/lib/')
from librgm.plotid import plotid


def plot_compass_arrow(AstWCS=None, direction='NS',
                       location=None,
                       length=2.0,
                       overplot=False):
    """

    Matplotlib arrow support:

    http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.arrow

    e.g.

    plt.arrow(x, y, dx, dy,
              head_width=0.05, head_length=0.1,
              facecolor='k', edgecolor='k')


    # Object oriented method
    # see https://scipy.github.io/old-wiki/pages/Cookbook/Matplotlib/Arrows.html
    Get the subplot that we are currently working on
    ax = gca()

    # Now add the arrow
    ax.add_patch(arr)


    http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.annotate

    http://matplotlib.org/api/patches_api.html#matplotlib.patches.ArrowStyle


    location code examples are here:

    http://matplotlib.org/api/legend_api.html#matplotlib.legend.Legend


    Astropy

    http://docs.astropy.org/en/stable/nddata/utils.html

    http://docs.astropy.org/en/stable/nddata/utils.html#d-cutout-from-a-skycoord-position


    """

    import numpy as np

    from matplotlib import pyplot as plt
    from matplotlib.patches import Circle, Ellipse

    demo = True

    print('Compute the on-sky position angle (East of North)')

    c1 = SkyCoord(0*u.deg, 0*u.deg)
    Separators = ' '
    precision = 1
    print('c1:',
          c1.ra.to_string(unit=u.hour,
          precision=precision+1, sep=' ', pad=True),
          c1.dec.to_string(unit=u.degree,
          precision=precision, sep=' ',
          pad=True, alwayssign=True))

    c2 = SkyCoord(1*u.deg, 0*u.deg)

    print('c2:',
          c2.ra.to_string(unit=u.hour,
          precision=precision+1, sep=' ', pad=True),
          c2.dec.to_string(unit=u.degree,
          precision=precision, sep=' ',
          pad=True, alwayssign=True))

    PA = c1.position_angle(c2).degree
    print('PA:', PA)

    c1 = SkyCoord(1*u.deg, 0*u.deg)
    c2 = SkyCoord(0*u.deg, 0*u.deg)

    print('c1:',
          c1.ra.to_string(unit=u.hour,
          precision=precision+1, sep=' ', pad=True),
          c1.dec.to_string(unit=u.degree,
          precision=precision, sep=' ',
          pad=True, alwayssign=True))

    print('c2:',
          c2.ra.to_string(unit=u.hour,
          precision=precision+1, sep=' ', pad=True),
          c2.dec.to_string(unit=u.degree,
          precision=precision, sep=' ',
          pad=True, alwayssign=True))


    # Computes the on-sky position angle (East of North)
    PA = c1.position_angle(c2).degree
    print('PA:', PA)

    if not overplot:
        plt.figure(figsize=(8,6))

        plt.plot([-10.0,10.0], [-10.0, 10.0], '.r')

    xrange = 20
    yrange = 20

    xarrow = 0.0
    yarrow = 0.0
    dxarrow = 0.0
    dyarrow = 2.0

    # plt.arrow( x, y, dx, dy, **kwargs )
    plt.arrow(xarrow, yarrow, dxarrow, dyarrow,
              head_width=0.4, head_length=0.4,
              facecolor='k', edgecolor='k')

    plt.arrow(0.0, 0.0, 2.0, 0.0,
              head_width=0.4, head_length=0.4,
              facecolor='k', edgecolor='k')

    plt.annotate('E', xy=(2.5, 0.0),
                 xycoords='data',
                 horizontalalignment='left',
                 verticalalignment='center',
                 fontsize='large')

    plt.annotate('N', xy=(0.0, 2.5),
                 xycoords='data',
                 horizontalalignment='center',
                 verticalalignment='bottom',
                 fontsize='large')

    plt.plot([-2.5, 2.5], [-8.0, -8.0], 'k', linewidth=2)
    plt.annotate('5"', xy=(0.0, -8.5),
                 xycoords='data',
                 horizontalalignment='center',
                 verticalalignment='top',
                 fontsize='large')

    plt.suptitle('On-sky position angle (PA) (East of North)')

    plt.xlabel('Arc seconds')
    plt.ylabel('Arc seconds')
    plt.axes().set_aspect('equal')

    # plot ellipse
    x0 = 0.0
    y0 = 0.0
    width = 4
    height = 8
    angle = 45.0
    theta = np.arange(0.0, 360.0, 1.0) * (np.pi / 180.0)
    xdata = x0 + (0.5 * width * np.cos(theta))
    ydata = y0 + (0.5 * height * np.sin(theta))

    rtheta = np.radians(angle)
    R = np.array([
                  [np.cos(rtheta), -np.sin(rtheta)],
                  [np.sin(rtheta),  np.cos(rtheta)],
                 ])

    xdata, ydata = np.dot(R, np.array([xdata, ydata]))

    plt.plot(xdata, ydata, 'b')

    plt.annotate('PA = ' + str(angle),
                 xy=(-9.0, 9.0),
                 xycoords='data',
                 horizontalalignment='left',
                 verticalalignment='center',
                 fontsize='large')

    xcenter = 5
    ycenter = 5
    width = 2
    height = width * 3
    angle = 0
    e1 = patches.Ellipse((xcenter, ycenter), width, height,
                         color = 'b',
                         angle=angle, linewidth=1,
                         fill=False, zorder=2)

    ax = plt.gca()
    ax.add_patch(e1)

    angle = 30
    e2 = patches.Ellipse((xcenter, ycenter), width, height,
                         color = 'r',
                         angle=angle, linewidth=1,
                         fill=False, zorder=2)

    ax = plt.gca()
    ax.add_patch(e2)


    angle = 60
    e2 = patches.Ellipse((xcenter, ycenter), width, height,
                         color = 'r',
                         angle=angle, linewidth=1,
                         fill=False, zorder=2)

    ax = plt.gca()
    ax.add_patch(e2)


    angle = 90
    e2 = patches.Ellipse((xcenter, ycenter), width, height,
                         color = 'r',
                         angle=angle, linewidth=1,
                         fill=False, zorder=2)

    ax = plt.gca()
    ax.add_patch(e2)


    plt.grid()
    plotid()
    plt.show()


def plot_scale_bar():
    """

    """



def radec2xieta(ra, dec, racent, deccent, xi, eta):
    """
    ra, dec in degrees

    Astrometry conversion from (alpha ,delta) to (xi ,eta)

    Standard coordinate conversion (see Smart, p283)

    Gnomonic projection or Tangent plane

    based on astrd2sn.pro

    """

    import numpy as np

    beta      = ra - racen
    cosbeta   = np.cos(np.deg2rad(beta))

    tandec    = np.tan(np.deg2rad(dec))
    tandeccen = np.tan(np.deg2rad(deccen))

    s   = cosbeta + tandec*tandeccen

    xi  = np.sin(np.deg2rad(beta))/np.cos(np.deg2rad(deccen))/s
    eta = (tandec - tandeccen*cosbeta)/s

    return xi, eta


def eq_to_tan(ra, dec, ra0=0.0, dec0=0.0):
    """Converts RA,Dec coordinates to xi, eta tangential coordiantes.
    See Olkin:1996 eq 3 for example, or Smart 1977.
    :return: tuple of xi, eta in degrees.

    """
    r = ra * np.pi / 180.
    d = dec * np.pi / 180.
    r0 = ra0 * np.pi / 180.
    d0 = dec0 * np.pi / 180.

    xi = np.cos(d) * np.sin(r - r0) \
        / (np.sin(d0) * np.sin(d)
        + np.cos(d0) * np.cos(d) * np.cos(r-r0))

    eta = (np.cos(d0) * np.sin(d)
            - np.sin(d0) * np.cos(d) * np.cos(r - r0)) \
    / (np.sin(d0) * np.sin(d) + np.cos(d0) * np.cos(d) * np.cos(r - r0))

    xi = xi * 180. / np.pi
    eta = eta * 180. / np.pi
    return xi, eta

def tan_to_eq(xiDeg, etaDeg, ra0Deg=0.0, dec0Deg=0.0):
    """
    Convert tangential coordinates to equatorial (RA, Dec) in degrees

    https://raw.githubusercontent.com/jonathansick/andromap/master/andromap/tanproj.py

    """

    xi = xiDeg * np.pi / 180.
    eta = etaDeg * np.pi / 180.
    ra0 = ra0Deg * np.pi / 180.
    dec0 = dec0Deg * np.pi / 180.

    ra = np.arctan(xi / (np.cos(dec0) - eta * np.sin(dec0))) + ra0
    dec = np.arctan((np.sin(dec0) + eta * np.cos(dec0))
            / (np.cos(dec0) - eta * np.sin(dec0))) * np.cos(ra - ra0)

    ra = ra * 180. / np.pi
    dec = dec * 180. / np.pi

    return ra, dec


def source_info(ra_source, dec_source,
    ra_table, dec_table, data=None, SearchRadius=4.0, format='COADD'):
    """
    uses astropy for the matching

    """

    print('Number of rows in table:', len(data))
    print('Format:', format)
    print('ra_source: ', ra_source)
    print('dec_source:', dec_source)
    ra_source = ra_source*u.degree
    dec_source = dec_source*u.degree
    print(Angle(ra_source), Angle(dec_source))

    # put source into list rather than a scalar
    radec_source = SkyCoord(Angle([ra_source]), Angle([dec_source]))
    print(radec_source)

    catalog=SkyCoord(ra=ra_table*u.degree,
        dec=dec_table*u.degree)
    print(catalog[0])

    idx, d2d, d3d = radec_source.match_to_catalog_sky(catalog)

    idxsource, idxcatalog, d2d, d3d = \
        catalog.search_around_sky(radec_source, SearchRadius*u.arcsec)

    print('Number of matched sources:', len(idxsource), len(idxcatalog))

    print(idxsource)

    print(idxcatalog)

    print(d2d.arcsec)

    imatch=-1
    if format == 'COADD':
        print('sep, MAG_PSF_I, MAG_AUTO_I')
        for index in idxcatalog:
            imatch=imatch+1
            print(
                data['COADD_OBJECTS_ID'][index],
                "{0:.3f}".format(d2d[imatch].arcsec),
                "{0:.2f}".format(data['MAG_PSF_I'][index]),
                "{0:.2f}".format(data['MAG_AUTO_I'][index]),
                data['OBJECT_NUMBER'][index])

    if format == 'SE':
        isort=np.argsort(data['MJD_OBS'][idxcatalog])
        print('Sort by MJD_OBS')
        print('NITE, MJD_OBS, BAND, CATALOGID, CCD, COADD_OBJECTS_ID, ' +
            'MAG_PSF, ZEROPOINT, CAL(MAG_PSF), FWHM')
        for index in isort:
            print(np.char.strip(data['NITE'][idxcatalog][index]),
                "{0:.7f}".format(data['MJD_OBS'][idxcatalog][index]),
                np.char.strip(data['BAND'][idxcatalog][index]),
                data['CATALOGID'][idxcatalog][index],
                data['CCD'][idxcatalog][index],
                data['COADD_OBJECTS_ID'][idxcatalog][index],
                "{0:.3f}".format(d2d[index].arcsec),
                "{0:.2f}".format(data['MAG_PSF'][idxcatalog][index]),
                "{0:.2f}".format(data['ZEROPOINT'][idxcatalog][index]),
                "{0:.2f}".format(data['MAG_PSF'][idxcatalog][index]),
                data['ZEROPOINT'][idxcatalog][index] +
                data['MAG_PSF'][idxcatalog][index],
                data['FWHM'][idxcatalog][index])


        print()
        for id in idxcatalog:
            imatch=imatch+1
            print(np.char.strip(data['NITE'][id]),
                "{0:.7f}".format(data['MJD_OBS'][id]),data['BAND'][id],
                "{0:.3f}".format(d2d[imatch].arcsec),
                data['CATALOGID'][id],
                data['COADD_OBJECTS_ID'][id],
                data['IMAGEPATH'][id],
                data['CATALOGPATH'][id])


def plot_lightcurve():
    """

    """


    return


def get_cutout(infile=None, data=None, ext=1, header=None,
               AstWCS=None,
               position=None, format='pixels', size=100,
               title=None, suptitle=None,
               imagetype='data',
               segmap=False, weightmap=False,
               plot=False, saveplot=True,
               plotfile_suffix=None, plotfile_prefix=None,
               verbose=False, debug=False):
    """



    """

    from astropy.nddata.utils import Cutout2D

    print('position: ', position[0], position[1])
    position = np.rint(position)
    print('position: ', position[0], position[1])

    if infile is not None:
        hdulist = fits.open(infile)
        hdulist.info()
        AstWCS = wcs.WCS(hdulist[ext].header)
        xpix0 = np.rint(position[0]) - (size/2)
        ypix0 = np.rint(position[1]) - (size/2)
        xpix0 =  xpix0.astype(int)
        ypix0 =  ypix0.astype(int)
        print('xpix0, ypix0: ', xpix0, ypix0)
        xpix1 = xpix0 + size
        ypix1 = ypix0 + size
        data = hdulist[ext].data[ypix0:ypix1, xpix0:xpix1]

    if debug: help(data)
    print('data.shape: ', data.shape)
    median=np.median(data)
    print('median.shape: ', median.shape)
    print('median: ', median)

    if segmap:
        # determine the list of unique sources in the segmentation image
        unique_sources = np.unique(data)
        nsources = len(unique_sources)
        print('Number of unique segmented sources: ', nsources)
        print(unique_sources)
        isource = 1
        # skip the first with value = zero which is background
        for unique_source in unique_sources[1:]:
            isource = isource + 1
            print(isource, unique_source)
            index = np.where(data == unique_source)
            print(index)
            data[index] = isource

    if ext != 2:
        itest = data > 0.5
        print('min: ', np.min(data[itest]))
        threshold = np.min(data[itest]) - 1

        print('threshold: ', threshold)

    print('max: ', np.max(data))

    mad  = median_absolute_deviation(data)
    mad_stdev = mad_std(data)

    print('mad:', mad)
    print('mad_std:', mad_stdev)
    print('mad/mad_std:', mad_stdev/mad)

    if ext != 2:
        data = data - threshold
        itest = data < 0
        data[itest] = 0

    median=np.median(data)
    print('median: ', median)

    position = (size/2, size/2)
    cutout = Cutout2D(data, position, size)
    if debug: help(cutout)

    if plot:

        plt.figure(figsize=(8,6))

        cmap = mpl.cm.jet
        if segmap:
            #cmap = mpl.cm.jet_r
            #cmap.set_under(color='w')
            #cmax = np.max(data)
            #cmap.set_clim(1,cmax)
            #itest = data > 0.5
            #data[itest] = np.nan
            data = ma.masked_where(data < 0.5, data)
            cmap.set_bad('w')
        #plt.imshow(cutout.data, origin='lower', interpolation='nearest')
            plt.imshow(data, origin='lower', interpolation='nearest',
                cmap=cmap)

        if not segmap:
            crange = 50
            if weightmap:crange = 10
            lower = -1.0
            vmin = median+(lower*mad_stdev)
            vmax=  min([median+(crange*mad_stdev),max])
            plt.imshow(data, origin='lower', interpolation='nearest',
                cmap=cmap,
                vmin=vmin, vmax=vmax)

        plt.gca().invert_xaxis()

        plt.xlabel('pixels')
        plt.ylabel('pixels')
        if title is not None: plt.title(title)
        if suptitle is not None: plt.suptitle(suptitle)
        plt.colorbar()
        plotid()


        if saveplot:
            plotfile = 'cutout'
            if segmap:
                plotfile = 'cutout_segmap'
            if weightmap:
                plotfile = 'cutout_weightmap'

            if plotfile_suffix is not None:
                plotfile = plotfile + '_' + plotfile_suffix

            if plotfile_prefix is not None:
                plotfile = plotfile_prefix + '_' + plotfile

            plotfile = plotfile + '.png'
            print('Saving :', plotfile)
            plt.savefig(plotfile)

        plt.show()

    return cutout.data


def wcs_PixelScale(AstWCS, x, y, debug=False):
    """
    simple determination of the pixel scale by computing the change in
    RA, Dec by one pixel centered on a specific pixel

         4
      1  0  2
         3

    requires WCS object from astropy.wcs

    """

    ra1, dec1 = AstWCS.wcs_pix2world(x-0.5, y, 1)
    ra2, dec2 = AstWCS.wcs_pix2world(x+0.5, y, 1)

    if debug: print('dec1, dec2: ', dec1, dec2)
    dec = (dec1 + dec2)/2.0
    if debug: print('Dec, cos(Dec): ', dec, np.cos(np.deg2rad(dec)))
    RAScale = (ra1 - ra2) *3600 * np.cos(np.deg2rad(dec))
    if debug: print('RAScale:', RAScale)

    ra3, dec3 = AstWCS.wcs_pix2world(x, y-0.5, 1)
    ra4, dec4 = AstWCS.wcs_pix2world(x, y+0.5, 1)

    DecScale = (dec4 - dec3) *3600

    print('x, y, RAScale, DecScale:',
        x, y, RAScale, DecScale, RAScale/DecScale)

    return RAScale, DecScale


def get_wcs(infile, ext=1, verbose=False, debug=False):
    """
    read WSC from a FITs image and do some simple World <-> Pixel operations

    """

    if debug: help(wcs)

    # Load the FITS hdulist using astropy.io.fits
    hdulist = fits.open(infile)
    print(hdulist.info())

    # Parse the WCS keywords in the primary HDU
    #key=raw_input("Enter any key to continue: ")
    AstWCS = wcs.WCS(hdulist[ext].header)
    NAXIS1= hdulist[ext].header['NAXIS1']
    NAXIS2= hdulist[ext].header['NAXIS2']
    print('NAXIS1: ', NAXIS1)
    print('NAXIS2: ', NAXIS2)

    RAScale, DecScale = wcs_PixelScale(AstWCS, 5000.5, 5000.5)

    if debug: AstWCS.wcs.print_contents()

    #help(AstWCS)
    x1 = 1
    y1 = 1
    radec = AstWCS.wcs_pix2world(x1, y1, 1)
    print('RA1, Dec1:', radec)

    y2 = NAXIS1
    radec = AstWCS.wcs_pix2world(x1, y2, 1)
    print('RA1, Dec1:', radec)

    wcs_PixelScale(AstWCS, 1, 1)
    wcs_PixelScale(AstWCS, 1, NAXIS2)
    wcs_PixelScale(AstWCS, NAXIS1, 1)
    wcs_PixelScale(AstWCS, NAXIS1, NAXIS2)
    # compute pixel scale

    RAScale, DecScale = wcs_PixelScale(AstWCS, 5000.5, 5000.5)

    ra4, dec4 = AstWCS.wcs_pix2world(NAXIS1, NAXIS2, 1)
    print('RA4, Dec4: ', ra4, dec4)

    if debug: help(AstWCS)
    hdulist.close()

    footprint = AstWCS.calc_footprint()
    print('footprint: ', footprint)
    key=raw_input("Enter any key to continue: ")

    if verbose or debug:
        print(repr(hdulist[1].header))

        if debug: AstWCS.wcs.print_contents()
        if debug: help(AstWCS.wcs)

        print('NAME: ', AstWCS.wcs.name)
        print('NAXIS: ', AstWCS.wcs.naxis)
        #print('NAXIS1: ', AstWCS.wcs.naxes)
        #print('NAXIS1: ', AstWCS.wcs.naxis1)
        #print('NAXIS2: ', AstWCS.wcs.naxis2)
        print('CRPIX: ', AstWCS.wcs.crpix)
        print('CDELT: ', AstWCS.wcs.cdelt)
        print('CD: ', AstWCS.wcs.cd)
        print('CRVAL: ', AstWCS.wcs.crval)
        print('CTYPE: ', AstWCS.wcs.ctype)

        crpix1, crpix2 = AstWCS.wcs.crpix
        crval1, crval2 = AstWCS.wcs.crval

        print('CRPIX1, CRPIX2: ', crpix1, crpix2)
        ra, dec = AstWCS.wcs_pix2world(crpix1, crpix2, 1)
        print('RA, Dec: ', ra, dec)

        print('CRVAL1, CRVAL2: ', crval1, crval2)
        radec = ICRS(crval1*u.degree, crval2*u.degree)
        #print(radec.ra.to_string(u.degree))
        print('RA, Dec: ',
            radec.ra.to_string(unit=u.hour, sep=' ', precision=3),
            radec.dec.to_string(unit=u.degree, sep=' ', precision=2))

        xpix, ypix = AstWCS.wcs_world2pix(crval1, crval2, 1)
        print('X0, Y0: ', xpix, ypix)

    return AstWCS


def list_lenses():
    """

    """

    lenses = config.get('sources','lenses')
    print('lenses: ', lenses, len(lenses))
    lenses = lenses.split(',')
    print('lenses: ', lenses, len(lenses))
    ilens = -1
    for lens in lenses:
        ilens = ilens + 1
        print('lens: ', lens, len(lens))
        ra = config.get(lens,'ra')
        dec = config.get(lens,'dec')
        print('ilens, ra, dec: ', ilens, ra, dec)
        radec_format = config.get(lens,'radec_format')
        print('radec format:', radec_format)


    for lens in lenses:
        print(lens)



def edgedetection(image=None):
    """

    http://scikit-image.org/docs/dev/api/skimage.feature.html#skimage.feature.canny
    http://stackoverflow.com/questions/29434533/edge-detection-for-image-stored-in-matrix

    """

    import numpy as np
    from matplotlib import pyplot as plt
    from scipy import ndimage

    from skimage import feature
    from skimage.filters import roberts, sobel, scharr, prewitt

    image_save = image

    # for convenience
    im = image


    # Compute the Canny filter for two values of sigma
    # http://scikit-image.org/docs/dev/api/skimage.feature.html#skimage.feature.cannya

    print('Image pixel value range:', np.min(im), np.max(im))
    edges1 = feature.canny(im, sigma=0.0,
                           low_threshold=0.0, high_threshold=5.0)

    # display results
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(16, 6),
                                        sharex=True, sharey=True)

    ax1.imshow(im, cmap=plt.cm.gray)
    ax1.axis('off')
    ax1.set_title('noisy image', fontsize=20)

    ax2.imshow(edges1, cmap=plt.cm.gray)
    ax2.axis('off')
    ax2.set_title('Canny filter \n $\sigma=0.0$', fontsize=20)

    fig.tight_layout()

    plt.show()

    # could try sobel filter
    # http://scikit-image.org/docs/dev/auto_examples/plot_edge_filter.html

    image = im
    edge_roberts = roberts(image)
    edge_sobel = sobel(image)

    fig, (ax0, ax1) = plt.subplots(ncols=2,
                                   sharex=True, sharey=True,
                                   subplot_kw={'adjustable':'box-forced'})

    ax0.imshow(edge_roberts, cmap=plt.cm.gray)
    ax0.set_title('Roberts Edge Detection')
    ax0.axis('off')

    ax1.imshow(edge_sobel, cmap=plt.cm.gray)
    ax1.set_title('Sobel Edge Detection')
    ax1.axis('off')

    plt.tight_layout()

    plt.show()


    # http://www.scipy-lectures.org/advanced/image_processing/auto_examples/plot_find_edges.html
    """
    Finding edges with Sobel filters
    ==================================

    The Sobel filter is one of the simplest way of finding edges.
    """

    import numpy as np
    from scipy import ndimage
    import matplotlib.pyplot as plt

    im = np.zeros((256, 256))
    im[64:-64, 64:-64] = 1

    im = ndimage.rotate(im, 15, mode='constant')
    im = ndimage.gaussian_filter(im, 8)

    sx = ndimage.sobel(im, axis=0, mode='constant')
    sy = ndimage.sobel(im, axis=1, mode='constant')
    sob = np.hypot(sx, sy)

    plt.figure(figsize=(16, 6))
    plt.subplot(141)
    plt.imshow(im, cmap=plt.cm.gray)
    plt.axis('off')
    plt.title('square', fontsize=20)
    plt.subplot(142)
    plt.imshow(sx)
    plt.axis('off')
    plt.title('Sobel (x direction)', fontsize=20)
    plt.subplot(143)
    plt.imshow(sob)
    plt.axis('off')
    plt.title('Sobel filter', fontsize=20)

    im += 0.07*np.random.random(im.shape)

    sx = ndimage.sobel(im, axis=0, mode='constant')
    sy = ndimage.sobel(im, axis=1, mode='constant')
    sob = np.hypot(sx, sy)

    plt.subplot(144)
    plt.imshow(sob)
    plt.axis('off')
    plt.title('Sobel for noisy image', fontsize=20)

    plt.subplots_adjust(wspace=0.02, hspace=0.02,
                        top=1, bottom=0, left=0, right=0.9)

    plt.show()


    # now apply to our image
    im = image_save

    sx = ndimage.sobel(im, axis=0, mode='constant')
    sy = ndimage.sobel(im, axis=1, mode='constant')
    sob = np.hypot(sx, sy)

    plt.figure(figsize=(16, 6))
    plt.subplot(141)
    plt.imshow(im, cmap=plt.cm.gray)
    plt.axis('off')
    plt.title('square', fontsize=20)
    plt.subplot(142)
    plt.imshow(sx)
    plt.axis('off')
    plt.title('Sobel (x direction)', fontsize=20)
    plt.subplot(143)
    plt.imshow(sob)
    plt.axis('off')
    plt.title('Sobel filter', fontsize=20)

    sx = ndimage.sobel(im, axis=0, mode='constant')
    sy = ndimage.sobel(im, axis=1, mode='constant')
    sob = np.hypot(sx, sy)

    plt.subplot(144)
    plt.imshow(sob)
    plt.axis('off')
    plt.title('Sobel for noisy image', fontsize=20)

    plt.subplots_adjust(wspace=0.02, hspace=0.02,
                        top=1, bottom=0, left=0, right=0.9)

    plt.show()



    # http://www.scipy-lectures.org/advanced/image_processing/auto_examples/plot_clean_morpho.html
    # plot_clean_morpho.py


    # def make_image():
    # """
    #
    # """

    np.random.seed(1)
    n = 10
    l = 256
    im = np.zeros((l, l))
    points = l*np.random.random((2, n**2))
    im[(points[0]).astype(np.int), (points[1]).astype(np.int)] = 1
    im = ndimage.gaussian_filter(im, sigma=l/(4.*n))

    mask = (im > im.mean()).astype(np.float)


    img = mask + 0.3*np.random.randn(*mask.shape)

    binary_img = img > 0.5

    # Remove small white regions
    open_img = ndimage.binary_opening(binary_img)
    # Remove small black hole
    close_img = ndimage.binary_closing(open_img)

    plt.figure(figsize=(16, 4))

    l = 128

    plt.subplot(141)
    plt.imshow(binary_img[:l, :l], cmap=plt.cm.gray)
    plt.axis('off')
    plt.subplot(142)
    plt.imshow(open_img[:l, :l], cmap=plt.cm.gray)
    plt.axis('off')
    plt.subplot(143)
    plt.imshow(close_img[:l, :l], cmap=plt.cm.gray)
    plt.axis('off')
    plt.subplot(144)
    plt.imshow(mask[:l, :l], cmap=plt.cm.gray)

    plt.contour(close_img[:l, :l], [0.5], linewidths=2, colors='r')

    plt.axis('off')

    plt.subplots_adjust(wspace=0.02, hspace=0.3,
                        top=1, bottom=0.1, left=0, right=1)

    plt.show()

    # http://stackoverflow.com/questions/1560424/how-can-i-get-the-x-y-values-of-the-line-that-is-ploted-by-a-contour-plot-mat
    # http://stackoverflow.com/questions/5666056/matplotlib-extracting-data-from-contour-lines


     # from http://stackoverflow.com/questions/29434533/edge-detection-for-image-stored-in-matrix
    thresh1 = 1
    thresh2 = 2

    #Load image
    # im = sp.misc.imread('jBD9j.png')

    im = image_save
    #Get threashold mask for different regions
    gryim = np.mean(im[:,:,0:2],2)
    region1 =  (thresh1 < gryim)
    region2 =  (thresh2 < gryim)
    nregion1 = ~ region1
    nregion2 = ~ region2

    #Plot figure and two regions
    fig, axs = plt.subplots(2,2)
    axs[0,0].imshow(im)
    axs[0,1].imshow(region1)
    axs[1,0].imshow(region2)

    #Clean up any holes, etc (not needed for simple figures here)
    #region1 = sp.ndimage.morphology.binary_closing(region1)
    #region1 = sp.ndimage.morphology.binary_fill_holes(region1)
    #region1.astype('bool')
    #region2 = sp.ndimage.morphology.binary_closing(region2)
    #region2 = sp.ndimage.morphology.binary_fill_holes(region2)
    #region2.astype('bool')

    #Get location of edge by comparing array to it's
    #inverse shifted by a few pixels
    shift = -2
    edgex1 = (region1 ^ np.roll(nregion1,shift=shift,axis=0))
    edgey1 = (region1 ^ np.roll(nregion1,shift=shift,axis=1))
    edgex2 = (region2 ^ np.roll(nregion2,shift=shift,axis=0))
    edgey2 = (region2 ^ np.roll(nregion2,shift=shift,axis=1))

    #Plot location of edge over image
    axs[1,1].imshow(im)
    axs[1,1].contour(edgex1,2,colors='r',lw=2.)
    axs[1,1].contour(edgey1,2,colors='r',lw=2.)
    axs[1,1].contour(edgex2,2,colors='g',lw=2.)
    axs[1,1].contour(edgey2,2,colors='g',lw=2.)




def explore_wise2des(data=None):
    """

    """


if __name__ == "__main__":


    # plot_compass_arrow(direction='NS', length=2.0)
    # key=raw_input("Enter any key to continue: ")


    import argparse

    import configparser

    from matplotlib import pyplot as plt

    config = configparser.RawConfigParser()
    config.read('VDESJ2325-5229_analysis.cfg')

    description = ''
    epilog = ''
    parser =  argparse.ArgumentParser(
        description=description, epilog=epilog,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    format_default = 'COADD'
    parser.add_argument("--format", default=format_default,
        dest='format', help="format as SE or COADD [Default]")

    release_default = 'Y1A1'
    parser.add_argument("--release", default=release_default,
        dest='release', help="release  as a string e.g. Y1A1 [Default]")

    parser.set_defaults(source='VDESJ2325-5229')
    parser.add_argument("--source", dest='source',
        help="source  as a string e.g. VDESJ2325-5229 [Default]")

    parser.set_defaults(band='i')
    parser.add_argument("--band",
        dest='band', help="waveband g, r, i, z, Y e.g. i [Default]")

    parser.set_defaults(cutout=False)
    parser.add_argument("--cutout", action='store_true',
        dest='cutout', help="cutout option")

    parser.set_defaults(debug=False)
    parser.add_argument("--debug", action='store_true',
        dest='debug', help="debug option")

    parser.set_defaults(list=False)
    parser.add_argument("--list", action='store_true',
        dest='list', help="list the lenses")

    parser.add_argument("--plotimages", action='store_true',
        default=False, help="plot images option")

    parser.add_argument("--showplots", action='store_true', default=False,
                        help="showplots option")

    # 38 pixels is 10 aecsec in DES Y1A1 image with
    # pixel size 0.267 arcsec/pixel
    parser.add_argument("--size_pixels", type=int, default=38,
                        help="Size in pixels")

    parser.add_argument("--edgedetection", action='store_true',
        default=False, help="Segmentation map edge detection  option")

    parser.set_defaults(verbose=False)
    parser.add_argument("--verbose", action='store_true',
        dest='verbose', help="verbose option")

    parser.set_defaults(verbose=False)
    parser.add_argument("--zoom", action='store_true',
                        default=False, help="zoom option")


    parser.set_defaults(xkcd=False)
    parser.add_argument("--xkcd", action='store_true',
        dest='xkcd', help="xkcd cartoon plot style")

    args = parser.parse_args()

    cutout = args.cutout

    BAND = args.band
    print('band: ', BAND)

    debug = args.debug
    print('debug: ', debug)

    format = args.format
    print('format: ', format)

    plotimages = args.plotimages

    release = args.release
    print('release: ', release)

    showplots = args.showplots
    print('showplots: ', showplots)

    source = args.source
    print('source: ', source)

    verbose = args.verbose
    print('verbose: ', verbose)

    zoom = args.zoom
    print('zoom: ', zoom)

    if args.xkcd: plt.xkcd()

    if args.list:
        list_lenses()
        sys.exit()

    # debug the list of lenses in the config file
    list_lenses()

    key=raw_input("Enter any key to continue: ")

    # cutout size in arc seconds
    size_cutout = 30.0

    # zoom size for ra, dec plots in arc seconds
    size_zoom = 6.0
    size_zoom = 10.0
    #size_zoom = 15.0

    size = 100
    size = 38
    size= args.size_pixels

    xrange=[-size/2.0, size/2.0]
    yrange=[-size/2.0, size/2.0]

    # colours to use for grizY plots
    colors=['blue','green','orange','red', 'maroon']

    RELEASE = 'Y1A1'

    print('RELEASE: ', RELEASE)
    print('SOURCE: ', source)

    datapath_desroot = config.get(source, 'datapath_desroot')

    if format.upper() == 'WISE2DES':
        datapath_wise2des = config.get(source, 'datapath_wise2des')


    datapath_cats = config.get(source, 'datapath_cats')
    # datapath_cats =   /data/des/VDESJ2325-5229/
    # datapath = '/home/rgm/soft/des/easyaccess/'

    radec_format = config.get(source, 'radec_format')

    TILENAME = config.get(source, 'tilename')

    ra = config.get(source,'ra')
    dec = config.get(source,'dec')
    radec_format = config.get(source, 'radec_format')
    print('ra, dec: ', ra, dec)
    print('radec format:', radec_format)

    if radec_format == 'hmsdms':
        coord = SkyCoord(ra + ' ' + dec, unit=(u.hourangle, u.degree))
        ra_source = coord.ra.degree
        dec_source = coord.dec.degree
        print(ra_source, dec_source)

        ra = ra_source
        dec = dec_source

    ra = float(ra)
    dec = float(dec)

    ra_source = ra
    dec_source = dec

    ra0 = ra
    dec0 = dec

    datapath = datapath_desroot + '/' + source + '/'

    filename_COADD = RELEASE + '_COADD_OBJECTS_' + source + '.fits'

    filename_SingleEpoch = config.get(source,'filename_se')
    # filename_SingleEpoch ='Y1A1_SE_' + source + '.fits'

    # ra =   ra + (1.0/3600.0)
    # dec =  dec - (1.0/3600.0)

    # from DES database
    # filename_SingleEpoch = 'Y1A1_FINALCUT_VDES2325-5229.fits'
    # filename_SingleEpoch = 'VDESJ2325-5229_Y1A1_SingleEpochFinalCut.fits'

    print('filename_SingleEpoch:', filename_SingleEpoch)
    print('TILENAME: ', TILENAME)

    # create rectanglar limits for box +-30" and the SQL fragment
    # use the WISE position for now

    # convert height in arcsecsonds to decimal degrees and define RA limits
    # for central dec which is OK for small offsets
    # add a litte for margin
    radec_size = 60.0
    # add a 1" margin
    radec_size = radec_size + 1.0

    print('radec_size: ', radec_size)

    delta_dec = radec_size/ 3600.0
    dec_min = dec - (delta_dec*0.5)
    dec_max = dec + (delta_dec*0.5)

    delta_ra = radec_size/3600.0
    print('Delta RA: ', delta_ra * 3600)
    delta_ra = delta_ra / np.cos(np.deg2rad(dec))
    print('Delta RA: ', delta_ra * 3600.0)
    ra_min = ra - (delta_ra*0.5)
    ra_max = ra + (delta_ra*0.5)

    print('RA, Dec: ', ra, dec)
    print('RA limits: ', ra_min, ra_max,
        (ra_max - ra_min)*3600.0,
        (ra_max - ra_min)*3600.0* np.cos(np.deg2rad(dec)))
    print('Dec limits: ', dec_min, dec_max, (dec_max - dec_min)*3600.0)

    # creat SQL fragment to get objects around the position
    print("WHERE \ \n " + \
          "    (ra BETWEEN " + "{0:.5f}".format(ra_min) + \
          " AND " + "{0:.5f}".format(ra_max) + ") \ \n" + \
          "    AND (dec BETWEEN " + "{0:.5f}".format(dec_min) + \
          " AND " + "{0:.5f}".format(dec_max) + ")\n")


    if format == 'SE':
        infile = datapath + filename_SingleEpoch

    if format == 'COADD':
        infile = datapath + filename_COADD

    if format.upper() == 'WISE2DES':
        tilename = 'DES0406-5414'
        suffix_wise2des = '_WISEfp_DEScat_WISE_match'
        infile = datapath_wise2des + tilename + '/' + \
            tilename + suffix_wise2des + '.fits'
        catfile = infile

    print('infile: ', infile)

    inpath = '/data/desardata/' + RELEASE + '/' + TILENAME + '/'
    print('inpath: ', inpath)


    print('args.edgedetection:', args.edgedetection)
    if args.edgedetection:

        ext = 1
        fzformat = True

        segmapname = TILENAME + '_' + BAND + '_seg.fits'
        if fzformat: segmapname = segmapname + '.fz'
        segmapfile = inpath + '/segmap/' + segmapname

        plotfile_suffix = BAND + '_' + str(size)
        plotfile_prefix = source
        suptitle = segmapname
        title = source

        wcs_image=get_wcs(segmapfile, verbose=True, debug=debug)
        xpix, ypix = wcs_image.wcs_world2pix(ra_source, dec_source, 1)
        print('RA, Dec, X, Y: ', ra_source, dec_source, xpix, ypix)

        key = raw_input("Enter any key to continue to edge detector: ")

        segmap = get_cutout(infile=segmapfile, ext=ext,
            title=title, suptitle=suptitle,
            position=(xpix, ypix), format='pixels', segmap=True,
            size=size, plot=True, saveplot=True,
            plotfile_prefix=plotfile_prefix,
            plotfile_suffix=plotfile_suffix,
            verbose=False, debug=debug)

        print(len(segmap), segmap.shape, len(segmap.shape), segmap.size)
        edgedetection(im=segmap)


    if format.upper() != 'WISE2DES':
        catfile = inpath + TILENAME + '_' + BAND + '_cat.fits'

    print('Read catfile: ', catfile)
    catdata = Table.read(catfile)
    print(catdata.colnames)
    catdata.info('stats')

    if plotimages:
        # fz format
        fzformat=True
        imagename = TILENAME + '_' + BAND + '.fits'
        if fzformat: imagename = imagename + '.fz'
        imagefile = inpath + imagename

        segmapname = TILENAME + '_' + BAND + '_seg.fits'
        if fzformat: segmapname = segmapname + '.fz'
        segmapfile = inpath + '/segmap/' + segmapname

        print('Read WCS for: ', imagefile)
        if debug: help(get_wcs)
        wcs_image=get_wcs(imagefile, verbose=True, debug=debug)
        #help(wcs_image)
        #print(wcs_image.wcs.name)
        xpix, ypix = wcs_image.wcs_world2pix(ra_source, dec_source, 1)
        print('RA, Dec, X, Y: ', ra_source, dec_source, xpix, ypix)

        filename = imagename
        plotfile_prefix = source
        suptitle = filename
        title = source

        plotfile_suffix = BAND + '_' + str(size)
        get_cutout(infile = imagefile, ext=1,
            title=title, suptitle=suptitle,
            position=(xpix, ypix), format='pixels',
            size=size, plot=True, saveplot=True,
            plotfile_prefix=plotfile_prefix,
            plotfile_suffix=plotfile_suffix,
            verbose=False, debug=debug)

        # weight map
        ext=2
        get_cutout(infile = imagefile, ext=ext,
            title=title, suptitle=suptitle,
            position=(xpix, ypix), format='pixels', weightmap=True,
            size=size, plot=True, saveplot=True,
            plotfile_prefix=plotfile_prefix,
            plotfile_suffix=plotfile_suffix,
            verbose=False, debug=debug)


        ext=1
        segmap = get_cutout(infile = segmapfile, ext=ext,
            title=title, suptitle=suptitle,
            position=(xpix, ypix), format='pixels', segmap=True,
            size=size, plot=True, saveplot=True,
            plotfile_prefix=plotfile_prefix,
            plotfile_suffix=plotfile_suffix,
            verbose=False, debug=debug)


        key=raw_input("Enter any key to continue to edge detector: ")


        key=raw_input("Enter any key to continue: ")


    print('format: ', format)
    print('Reading file: ', infile)
    data = Table.read(infile)
    print(data.colnames)
    data.info('stats')

    WAVEBANDS = ['G','R','I','Z','Y']

    if zoom:
        xrange=[-size_zoom/2.0, size_zoom/2.0]
        yrange=[-size_zoom/2.0, size_zoom/2.0]

    if format == 'SE':

        ra  = data['RA']
        dec = data['DEC']

        print('Number of rows: ', len(ra))

        rawin = data['ALPHAWIN_J2000']
        decwin = data['DELTAWIN_J2000']

        ra = rawin
        dec = decwin

        # change waveband to upper case and strip blanks
        # not very elegant but shows the principle
        obsband = data['BAND']
        wavebands_unique = np.unique(obsband)
        print('Observations wavebands: ', wavebands_unique)

        wavebands_unique = map(str.upper,wavebands_unique)
        print(wavebands_unique)

        wavebands_unique = map(str.strip,wavebands_unique)
        print(wavebands_unique)

        obsband = map(str.upper, obsband)
        obsband = map(str.strip, obsband)

        ra_min = np.min(ra)
        ra_max = np.max(ra)
        print('RA range: ', ra_min, ra_max)
        delta_ra = (ra - ra0) * 3600.0 * np.cos(np.deg2rad(dec0))

        dec_min = np.min(dec)
        dec_max = np.max(dec)
        print('Dec range: ', dec_min, dec_max)
        delta_dec = (dec - dec0)*3600.0

        plt.figure(figsize=(6,6))

        xdata= delta_ra
        ydata= delta_dec

        print(np.min(delta_ra), np.max(delta_ra))
        print(np.min(delta_dec), np.max(delta_dec))

        itest = (xdata > xrange[0]) & (xdata < xrange[1]) & \
            (ydata > yrange[0]) & (ydata < yrange[1])

        xdata = xdata[itest]
        ydata = ydata[itest]

        ndata = len(xdata)
        print('Number of rows: ', len(ra))
        plt.plot(delta_ra, delta_dec, '.',
            label=str(ndata))

        if zoom:
            plt.xlim(xrange)
            plt.ylim(yrange)

        plt.title(infile, fontsize='medium')
        plt.grid(True)

        plt.xlabel('Delta RA (arc seconds)')
        plt.ylabel('Delta Dec (arc seconds)')

        plotid(progname=True)

        #plt.legend()
        #plt.show()


        iband = -1
        print('Cycle through wavebands')
        for WAVEBAND in WAVEBANDS:
            iband = iband + 1
            #print(WAVEBAND, obsband[0], len(obsband))
            #if iband ==1: help(obsband)
            # convert python list to numpy array
            obsband = np.asarray(obsband)
            itest = (obsband == WAVEBAND)
            print('itest: ', itest)
            print(WAVEBAND, len(obsband[itest]))
            print(obsband[itest])

            ra = data['ALPHAWIN_J2000'][itest]
            dec = data['DELTAWIN_J2000'][itest]

            print(len(ra))

            delta_ra = (ra - ra0)*3600.0 * np.cos(np.deg2rad(dec0))

            delta_dec = (dec - dec0)*3600.0

            xdata= delta_ra
            ydata= delta_dec
            print(np.min(delta_ra), np.max(delta_ra))
            print(np.min(delta_dec), np.max(delta_dec))

            itest = (xdata > xrange[0]) & (xdata < xrange[1]) & \
                (ydata > yrange[0]) & (ydata < yrange[1])

            xdata = xdata[itest]
            ydata = ydata[itest]

            print(np.min(xdata), np.max(xdata))
            print(np.min(ydata), np.max(ydata))

            ndata = len(xdata)

            print(iband, colors[iband], ndata)

            plt.plot(delta_ra, delta_dec, '.', color=colors[iband],
                label=WAVEBAND+': ' + str(ndata))

        plt.legend()

        plotfile = source + '_SingleEpoch_radec_zoom.png'
        plt.savefig(plotfile)
        #plt.clf()
        print('Saving: ', plotfile)

        plt.show()

        plt.close()

        # figsize units are w, h in inchs
        fig = plt.figure(figsize=(8, 8))
        #fig.subplots_adjust(bottom=0.025, left=0.025, top = 0.975, right=0.975)

        # loop through the wavebands and epochs
        ra  = data['RA']
        dec = data['DEC']

        delta_ra = (ra - ra0) * 3600.0 * np.cos(np.deg2rad(dec0))
        delta_dec = (dec - dec0) * 3600.0

        itest = (delta_ra > xrange[0]) & (delta_ra < xrange[1]) & \
            (delta_dec > yrange[0]) & (delta_dec < yrange[1])

        data = data [itest]
        isort =np.argsort(data['MJD_OBS'])
        CATALOGID_unique, indices = \
            np.unique(data['CATALOGID'], return_index=True)

        print('Number of unique CATALOGIDs: ', len(CATALOGID_unique))
        key=raw_input("\nEnter any key to continue: ")

        index = 0
        iplot = 0
        for CATALOGID in CATALOGID_unique:
            index = index + 1
            iplot = iplot + 1
            plt.subplot(3, 3, iplot)

            idata = (data['CATALOGID'] == CATALOGID)

            ra = data['ALPHAWIN_J2000'][idata]
            dec = data['DELTAWIN_J2000'][idata]
            WAVEBAND = data['BAND'][idata]

            print(len(ra))

            delta_ra = (ra - ra0)*3600.0 * np.cos(np.deg2rad(dec0))

            delta_dec = (dec - dec0)*3600.0

            print(np.min(delta_ra), np.max(delta_ra))
            print(np.min(delta_dec), np.max(delta_dec))


            itest = (delta_ra > xrange[0]) & (delta_ra < xrange[1]) & \
                (delta_dec > yrange[0]) & (delta_dec < yrange[1])

            xdata= delta_ra[itest]
            ydata= delta_dec[itest]
            WAVEBAND = WAVEBAND[itest][0].upper()
            WAVEBAND = WAVEBAND.strip()

            print(np.min(xdata), np.max(xdata))
            print(np.min(ydata), np.max(ydata))

            ndata = len(xdata)
            iband = WAVEBANDS.index(WAVEBAND)
            print(iband, colors[iband], ndata)

            plt.plot(delta_ra, delta_dec, '.', color=colors[iband],
                label=WAVEBAND+': ' + str(ndata))

            #plt.legend(numpoints=0)
            WAVEBAND=WAVEBAND.strip()
            title = str(CATALOGID) + ' ' + WAVEBAND
            plt.title(title, fontsize='small')
            plt.suptitle(source, fontsize='medium')

            #plt.text(0.5,0.5,'plt.text; Hello World',
            #    transform=plt.gca().transAxes)

            plt.annotate(WAVEBAND, (0.1, 0.8),
               textcoords='axes fraction')

            noticks=False
            if noticks:
                plt.xticks(())
                plt.yticks(())

            plt.grid(True)

            if zoom:
                plt.xlim(xrange)
                plt.ylim(yrange)


        plotfile = source + '_multi_SingleEpoch_radec_zoom.png'
        plt.savefig(plotfile)
        #plt.clf()
        print('Saving: ', plotfile)

        plt.show()
        plt.close()

        key=raw_input("\nEnter any key to continue: ")

    if format == 'COADD':

        itest=np.unique(data['TILENAME'])
        print(itest)
        ra = data['RA']
        dec = data['DEC']
        COADD_OBJECTS_ID = data['COADD_OBJECTS_ID']

        ra_min = np.min(ra)
        ra_max = np.max(ra)
        delta_ra = (ra - ra0)*3600.0 * np.cos(np.deg2rad(dec0))

        dec_min = np.min(dec)
        dec_max = np.max(dec)
        delta_dec = (dec - dec0)*3600.0

        xdata= delta_ra
        ydata= delta_dec

        print(xrange)
        print(yrange)
        itest = (xdata > xrange[0]) & (xdata < xrange[1]) & \
                (ydata > yrange[0]) & (ydata < yrange[1])

        xdata = xdata[itest]
        ydata = ydata[itest]
        ndata = len(xdata)

        plt.figure(figsize=(8,8))
        plt.axes().set_aspect('equal')
        plt.plot(delta_ra, delta_dec, '.k', label='COADD: '+str(ndata))

        if zoom:
            plt.xlim(xrange)
            plt.ylim(yrange)

        plt.title(infile, fontsize='medium')
        plt.grid(True)

        WAVEBANDS = ['G','R','I','Z','Y']

        print("id, object_id, WAVEBAND, NEPOCHS, MAG_AUTO, " + \
                  "Dra, Ddec, " + \
                  "FLUX_RADIUS, KRON_RADIUS, " + \
                  "A, B, PA, Aspect Ratio(A/B)")

        iband = -1
        for WAVEBAND in WAVEBANDS:
            iband = iband + 1
            ra = data['ALPHAWIN_J2000_'+WAVEBAND]
            dec = data['DELTAWIN_J2000_'+WAVEBAND]

            # used to convert pixels to arc seconds
            PIXEL_SIZE=0.27

            delta_ra = (ra - ra0)*3600.0 * np.cos(np.deg2rad(dec0))
            delta_dec = (dec - dec0)*3600.0

            xdata= delta_ra
            ydata= delta_dec

            # limit the data
            itest = (xdata > xrange[0]) & (xdata < xrange[1]) & \
                (ydata > yrange[0]) & (ydata < yrange[1])

            xdata = xdata[itest]
            ydata = ydata[itest]
            ndata = len(xdata)

            delta_ra = delta_ra[itest]
            delta_dec = delta_dec[itest]

            # maybe combine this with the assignment above
            COADD_OBJECTS_ID = data['COADD_OBJECTS_ID'][itest]
            OBJECT_NUMBER = data['OBJECT_NUMBER'][itest]
            RADEC_SOURCE_BAND = data['RADEC_SOURCE_BAND'][itest]

            # Detection image
            A_IMAGE = data['A_IMAGE'][itest]*PIXEL_SIZE
            B_IMAGE = data['B_IMAGE'][itest]*PIXEL_SIZE
            THETA_IMAGE = data['THETA_IMAGE'][itest]

            KRON_RADIUS = data['KRON_RADIUS'][itest]*PIXEL_SIZE

            XMIN_IMAGE=data['XMIN_IMAGE'][itest]
            XMAX_IMAGE=data['XMAX_IMAGE'][itest]
            YMIN_IMAGE=data['YMIN_IMAGE'][itest]
            YMAX_IMAGE=data['YMAX_IMAGE'][itest]

            # Measurement image
            ALPHAWIN_J2000 = data['ALPHAWIN_J2000_'+WAVEBAND][itest]
            DELTAWIN_J2000 = data['DELTAWIN_J2000_'+WAVEBAND][itest]
            XWIN_IMAGE = data['XWIN_IMAGE_'+WAVEBAND][itest]
            YWIN_IMAGE = data['YWIN_IMAGE_'+WAVEBAND][itest]

            NEPOCHS = data['NEPOCHS_'+WAVEBAND][itest]
            MAG_AUTO = data['MAG_AUTO_'+WAVEBAND][itest]
            FLUX_RADIUS = data['FLUX_RADIUS_'+WAVEBAND][itest]*PIXEL_SIZE

            AWIN_IMAGE = data['AWIN_IMAGE_'+WAVEBAND][itest]*PIXEL_SIZE
            BWIN_IMAGE = data['BWIN_IMAGE_'+WAVEBAND][itest]*PIXEL_SIZE
            THETAWIN_IMAGE = data['THETAWIN_IMAGE_'+WAVEBAND][itest]

            ELLIP1MODEL_WORLD = data['ELLIP1MODEL_WORLD_'+WAVEBAND][itest]
            ELLIP2MODEL_WORLD = data['ELLIP2MODEL_WORLD_'+WAVEBAND][itest]

            FWHM_WORLD = data['FWHM_WORLD_'+WAVEBAND][itest]

            ISOAREA_WORLD = data['ISOAREA_WORLD_'+WAVEBAND][itest]

            FLAGS = data['FLAGS_'+WAVEBAND][itest]

            alpha=0.1
            i =-1

            for id in xdata:
                i = i + 1

                circle = Circle([delta_ra[i], delta_dec[i]], 0.25,
                    edgecolor='none', facecolor=colors[iband], alpha=alpha)
                plt.gca().add_patch(circle)

                width = AWIN_IMAGE[i]
                height = BWIN_IMAGE[i]
                angle = THETAWIN_IMAGE[i] + 90.0
                coadd_objects_id = COADD_OBJECTS_ID[i]

                print(i, coadd_objects_id,
                    OBJECT_NUMBER[i],
                    WAVEBAND,
                    RADEC_SOURCE_BAND[i],
                    "{:4d}".format(NEPOCHS[i]),
                    "{:8.2f}".format(MAG_AUTO[i]),
                    "{:6.2f}".format(delta_ra[i]),
                    "{:6.2f}".format(delta_dec[i]),
                    "{:6.2f}".format(FLUX_RADIUS[i]),
                    "{:6.2f}".format(KRON_RADIUS[i]),
                    "{:6.2f}".format(width),
                    "{:6.2f}".format(height),
                    "{:7.1f}".format(angle),
                    "{:6.3f}".format(width/height))

                print(i, coadd_objects_id,
                    OBJECT_NUMBER[i],
                    WAVEBAND,
                    "{:8.2f}".format(A_IMAGE[i]),
                    "{:8.2f}".format(B_IMAGE[i]),
                    "{:7.1f}".format(THETA_IMAGE[i]))

                print(i, coadd_objects_id,
                    OBJECT_NUMBER[i],
                    WAVEBAND,
                    "{:7.1f}".format(XWIN_IMAGE[i]),
                    "{:7.1f}".format(YWIN_IMAGE[i]),
                    "{:7.1f}".format(XMIN_IMAGE[i]),
                    "{:7.1f}".format(XMAX_IMAGE[i]),
                    "{:7.1f}".format(YMIN_IMAGE[i]),
                    "{:7.1f}".format(YMAX_IMAGE[i]),
                    "{:5.1f}".format(XMAX_IMAGE[i]-XMIN_IMAGE[i]),
                    "{:5.1f}".format(YMAX_IMAGE[i]-YMIN_IMAGE[i]),
                    "{:7.1f}".format((XMIN_IMAGE[i]+XMAX_IMAGE[i])/2.0),
                    "{:7.1f}".format((YMIN_IMAGE[i]+YMAX_IMAGE[i])/2.0),
                    "{:4d}".format(FLAGS[i]))

                ellipse = Ellipse([delta_ra[i], delta_dec[i]],
                    width=width/2.0, height=height/2.0, angle=angle,
                    edgecolor='none', facecolor=colors[iband], alpha=alpha)
                plt.gca().add_patch(ellipse)

            plt.plot(delta_ra, delta_dec, '.', color=colors[iband],
                label=WAVEBAND + ': ' + str(ndata))

        plt.suptitle('Windowed positions')

        plt.xlabel('Delta RA (arc seconds)')
        plt.ylabel('Delta Dec (arc seconds)')
        plt.legend(fontsize='medium')
        plotid(progname=True)

        plotfile = source + '_COADD_radec_zoom.png'
        plt.savefig(plotfile)
        #plt.clf()
        print('Saving: ', plotfile)

        plt.show()

    if format.upper() == 'WISE2DES':

        ra = data['RA_CALC_I']
        dec = data['DEC_CALC_I']

        ra_min = np.min(ra)
        ra_max = np.max(ra)
        delta_ra = (ra - ra0)*3600.0 * np.cos(np.deg2rad(dec0))

        dec_min = np.min(dec)
        dec_max = np.max(dec)
        delta_dec = (dec - dec0)*3600.0

        xdata= delta_ra
        ydata= delta_dec

        itest = (xdata > xrange[0]) & (xdata < xrange[1]) & \
                (ydata > yrange[0]) & (ydata < yrange[1])

        xdata = xdata[itest]
        ydata = ydata[itest]
        ndata = len(xdata)

        plt.figure(figsize=(8,8))
        plt.axes().set_aspect('equal')
        plt.plot(delta_ra, delta_dec, '.k', label='COADD: '+str(ndata))

        if zoom:
            plt.xlim(xrange)
            plt.ylim(yrange)

        plt.title(infile, fontsize='medium')
        plt.grid(True)



    # match the objects to the nominal Lensed source position

    if format == 'SE':
        ra_table  = data['RA']
        dec_table = data['DEC']

    if format == 'COADD':
        ra_table  = data['RA']
        dec_table = data['DEC']

    SearchRadius=4.0
    source_info(ra_source, dec_source,
        ra_table, dec_table, data=data,
        SearchRadius=SearchRadius,
        format=format)
