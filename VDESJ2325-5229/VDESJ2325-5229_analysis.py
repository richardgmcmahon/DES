
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

    print('Number of rows in table: ', len(data))
    print('Format: ', format)
    print('ra_source: ', ra_source)
    print('dec_source: ', dec_source)
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

    print('Number of matched sources: ', len(idxsource), len(idxcatalog))

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


def get_cutout(infile=None, data=None, header=None, WCS=None, ext=1,
    position=None, format='pixels', size=100, title=None, 
    segmap=False, weightmap=False,
    plot=False, saveplot=True,
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
        WCS = wcs.WCS(hdulist[ext].header)
        xpix0 = np.rint(position[0]) - (size/2)
        ypix0 = np.rint(position[1]) - (size/2)
        xpix0 =  xpix0.astype(int)
        ypix0 =  ypix0.astype(int)
        print('xpix0, ypix0: ', xpix0, ypix0)
        xpix1 = xpix0 + 100
        ypix1 = ypix0 + 100
        data = hdulist[ext].data[ypix0:ypix1, xpix0:xpix1]

    if debug: help(data)
    print('data.shape: ', data.shape)
    median=np.median(data)
    print('median.shape: ', median.shape)
    print('median: ', median)

    if segmap:
        unique_sources = np.unique(data)
        nsources = len(unique_sources)
        print('Number of unique segmented sources: ', nsources)
        print(unique_sources)
        isource = 1
        # skip the first with value = zero which isbackground
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

    print('mad: ', mad)
    print('mad_std: ', mad_stdev)
    print('mad/mad_std: ', mad_stdev/mad)

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
        cmap = mpl.cm.jet
        if segmap: 
            #cmap = mpl.cm.jet_r
            cmap.set_bad('w')
            #cmap.set_under(color='w')
            #cmax = np.max(data)
            #cmap.set_clim(1,cmax)
            #itest = data > 0.5
            #data[itest] = np.nan
            data = ma.masked_where(data < 0.5, data)
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

        plt.xlabel('pixels')
        plt.ylabel('pixels')
        if title is not None: plt.title(title)
        plt.colorbar()


        if saveplot:
            plotfile='cutout.png'
            if segmap: plotfile='cutout_segmap.png'
            if weightmap: plotfile='cutout_weightmap.png'
            print('Saving :', plotfile)
            plt.savefig(plotfile)

        plt.show()

    return cutout


def plot_cutout():


    return


def wcs_PixelScale(WCS, x, y, debug=False):
    """
    simple determination of the pixel scale by computing the change in
    RA, Dec by one pixel centered on a specific pixel

         4
      1  0  2
         3

    requires WCS object from astropy.wcs

    """

    ra1, dec1 = WCS.wcs_pix2world(x-0.5, y, 1)
    ra2, dec2 = WCS.wcs_pix2world(x+0.5, y, 1)

    if debug: print('dec1, dec2: ', dec1, dec2)
    dec = (dec1 + dec2)/2.0
    if debug: print('Dec, cos(Dec): ', dec, np.cos(np.deg2rad(dec)))
    RAScale = (ra1 - ra2) *3600 * np.cos(np.deg2rad(dec))
    if debug: print('RAScale:', RAScale)

    ra3, dec3 = WCS.wcs_pix2world(x, y-0.5, 1)
    ra4, dec4 = WCS.wcs_pix2world(x, y+0.5, 1)

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
    #WCS = wcs.WCS(hdulist[0].header)
    #key=raw_input("Enter any key to continue: ")
    WCS = wcs.WCS(hdulist[ext].header)
    NAXIS1= hdulist[ext].header['NAXIS1']
    NAXIS2= hdulist[ext].header['NAXIS2']
    print('NAXIS1: ', NAXIS1)
    print('NAXIS2: ', NAXIS2)
    ra1, dec1 = WCS.wcs_pix2world(1, 1, 1)  
    print('RA1, Dec1: ', ra1, dec1)
    wcs_PixelScale(WCS, 1, 1)
    wcs_PixelScale(WCS, 1, NAXIS2)
    wcs_PixelScale(WCS, NAXIS1, 1)
    wcs_PixelScale(WCS, NAXIS1, NAXIS2)
    # compute pixel scale

    RAScale, DecScale = wcs_PixelScale(WCS, 5000.5, 5000.5)

    ra4, dec4 = WCS.wcs_pix2world(NAXIS1, NAXIS2, 1)  
    print('RA4, Dec4: ', ra4, dec4)

    if debug: help(WCS)
    hdulist.close()

    footprint = WCS.calc_footprint()
    print('footprint: ', footprint)
    key=raw_input("Enter any key to continue: ")

    if verbose or debug:
        print(repr(hdulist[1].header))

        if debug: WCS.wcs.print_contents()
        if debug: help(WCS.wcs)

        print('NAME: ', WCS.wcs.name)
        print('NAXIS: ', WCS.wcs.naxis)
        #print('NAXIS1: ', WCS.wcs.naxes)
        #print('NAXIS1: ', WCS.wcs.naxis1)
        #print('NAXIS2: ', WCS.wcs.naxis2)
        print('CRPIX: ', WCS.wcs.crpix)
        print('CDELT: ', WCS.wcs.cdelt)
        print('CD: ', WCS.wcs.cd)
        print('CRVAL: ', WCS.wcs.crval)
        print('CTYPE: ', WCS.wcs.ctype)

        crpix1, crpix2 = WCS.wcs.crpix
        crval1, crval2 = WCS.wcs.crval

        print('CRPIX1, CRPIX2: ', crpix1, crpix2)
        ra, dec = WCS.wcs_pix2world(crpix1, crpix2, 1)  
        print('RA, Dec: ', ra, dec)

        print('CRVAL1, CRVAL2: ', crval1, crval2)
        radec = ICRS(crval1*u.degree, crval2*u.degree)
        #print(radec.ra.to_string(u.degree))
        print('RA, Dec: ',
            radec.ra.to_string(unit=u.hour, sep=' ', precision=3),
            radec.dec.to_string(unit=u.degree, sep=' ', precision=2))

        xpix, ypix = WCS.wcs_world2pix(crval1, crval2, 1)
        print('X0, Y0: ', xpix, ypix)

    return WCS


def list_lenses():
    """

    """

    lenses = ['VDESJ2325-5229','DESJ0408-5354',
        'DESJ0115-5244', 'DESJ2146-0047']
   
    """
    
    DES J0115-5244 01 15 57.32 -52:44:23.2
    
    DES J2146-0047 21 46 46.04 -00:47:44.3

    """

    for lens in lenses:
        print(lens)


if __name__ == "__main__":

    from argparse import ArgumentParser

    parser = ArgumentParser()

    format_default = 'COADD'
    parser.add_argument("--format", default=format_default,
        dest='format', help="format as SE or COADD [Default]")

    release_default = 'Y1A1'
    parser.add_argument("--release", default=release_default,
        dest='release', help="release  as a string e.g. Y1A1 [Default]")

    source_default = 'VDESJ2325-5229'
    parser.add_argument("--source", default=source_default,
        dest='source', 
        help="source  as a string e.g. DESJ0408-5354, VDESJ2325-5229 [Default]")

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

    parser.set_defaults(showplots=False)
    parser.add_argument("--showplots", action='store_true',
        dest='showplots', help="showplots option")

    parser.set_defaults(verbose=False)
    parser.add_argument("--verbose", action='store_true',
        dest='verbose', help="verbose option")

    parser.set_defaults(xkcd=False)
    parser.add_argument("--xkcd", action='store_true',
        dest='xkcd', help="xkcd cartoon plot style")

    args = parser.parse_args()

    cutout = args.cutout

    debug = args.debug
    print('debug: ', debug)

    format = args.format
    print('format: ', format)

    release = args.release
    print('release: ', release)

    showplots = args.showplots
    print('showplots: ', showplots)

    source = args.source
    print('source: ', source)

    verbose = args.verbose
    print('verbose: ', verbose)

    if args.xkcd: plt.xkcd()

    if args.list: 
        list_lenses()
        sys.exit()

    # cutout size in arc seconds
    cutout=True
    size_cutout = 30.0


    zoom = True
    # zoom size for ra, dec plots in arc seconds
    size_zoom = 6.0
    #size_zoom = 15.0

    # colours to use for grizY plots
    colors=['blue','green','orange','red', 'maroon']

    RELEASE = 'Y1A1'

    print('RELEASE: ', RELEASE)
    print('SOURCE: ', source)

    if source == 'DESJ0115-5244':

        TILENAME = 'DES0115-5248'

        ra='01 15 57.32'
        dec='-52 44 23.2'

        coord = SkyCoord(ra+' '+dec, unit=(u.hourangle, u.degree))
        ra_source = coord.ra.degree
        dec_source = coord.dec.degree

        print(ra_source, dec_source)
  
        ra = ra_source
        dec = dec_source

        datapath = '/data/des/' + source + '/'

        filename_COADD = RELEASE + '_COADD_OBJECTS_' + source + '.fits'

    if source == 'DESJ0408-5354':

        TILENAME='DES0406-5414'

        ra_wise =   62.090135 
        dec_wise = -53.89970 

        ra= '04h08m21.9199s' 
        dec= '-53d54m00.9576s'

        ra =   62.090135 + (1.0/3600.0)
        dec = -53.89970 - (1.0/3600.0)

        ra_source=ra_wise
        dec_source=dec_wise

        ra0 = ra
        dec0 = dec

        datapath = '/home/rgm/soft/des/easyaccess/'

        filename_SingleEpoch ='Y1A1_SE_hlin.fits'

        filename_COADD = RELEASE + '_COADD_OBJECTS_' + source + '.fits'

    if source == 'DESJ2146-0047':

        TILENAME = 'DES2145-0041'

        ra = '21 46 46.0' 
        dec = '-00 47 44.3'

        coord = SkyCoord(ra+' '+dec, unit=(u.hourangle, u.degree))
        ra_source = coord.ra.degree
        dec_source = coord.dec.degree

        print(ra_source, dec_source)
  
        ra = ra_source
        dec = dec_source

        datapath = '/data/des/' + source + '/'

        filename_COADD = 'Y1A1_COADD_OBJECTS_' + source + '.fits'

    if source == 'VDESJ2325-5229':

        TILENAME='DES2327-5248'

        #ra_wise = '351.4217576 degrees'	
        #dec_wise = '-52.4875349 degrees'

        ra_wise = 351.4217576
        dec_wise = -52.4875349

        ra='23:25:41.20 hours' 
        dec='-52:29:15.2 degrees'

        ra = 351.4217576
        dec = -52.4875349

        ra0 = 351.4217576
        dec0 = -52.4875349

        ra_source = ra_wise
        dec_source = dec_wise


        datapath = '/data/des/VDES2325-5229/'

        # from DES database
        filename_SingleEpoch = 'Y1A1_FINALCUT_VDES2325-5229.fits'
        filename_SingleEpoch = 'VDESJ2325-5229_Y1A1_SingleEpochFinalCut.fits'

        filename_COADD = 'Y1A1_COADD_OBJECTS_VDES2325-5229.fits'

    
    print('TILENAME: ', TILENAME)




    # create rectanglar limits for box +-30" and the SQL fragment
    # use the WISE position for now

    # convert height in arcsecsonds to decimal degrees and define RA limits
    # for central dec which is OK for small offsets
    # add a litte for margin
    radec_size = 60.0 + 1.0
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


    inpath = '/data/desardata/' + RELEASE + '/' + TILENAME + '/'

    BAND = 'i'
    catfile = inpath + TILENAME + '_' + BAND + '_cat.fits'
    print('Read catfile: ', catfile)
    catdata = Table.read(catfile)
    print(catdata.colnames)
    catdata.info('stats')

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
    get_cutout(infile = imagefile, ext=1, title=filename,
        position=(xpix, ypix), format='pixels',
        size=100, plot=True, saveplot=True,
        verbose=False, debug=debug)


    # weight map
    ext=2
    get_cutout(infile = imagefile, ext=ext, title=filename,
        position=(xpix, ypix), format='pixels', weightmap=True,
        size=100, plot=True, saveplot=True,
        verbose=False, debug=debug)


    filename = segmapname
    ext=1
    get_cutout(infile = segmapfile, ext=ext, title=filename,
        position=(xpix, ypix), format='pixels', segmap=True,
        size=100, plot=True, saveplot=True,
        verbose=False, debug=debug)

    key=raw_input("Enter any key to continue: ")




    print('format: ', format)
    print('Reading file: ', infile)
    data = Table.read(infile)
    print(data.colnames)
    data.info('stats')


    WAVEBANDS = ['G','R','I','Z','Y']

    zoom = True
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
        delta_ra = (ra - ra0)*3600.0 * np.cos(np.deg2rad(dec0))

        dec_min = np.min(dec)
        dec_max = np.max(dec)
        print('Dec range: ', dec_min, dec_max)
        delta_dec = (dec - dec0)*3600.0

        plt.figure(figsize=(8,8))

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


        plt.xlabel('Delta RA (arc seconds)')
        plt.ylabel('Delta Dec (arc seconds)')
        plt.legend(fontsize='medium')
        plotid(progname=True)

        plotfile = source + '_COADD_radec_zoom.png'
        plt.savefig(plotfile)
        #plt.clf()
        print('Saving: ', plotfile)

        plt.show()

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


