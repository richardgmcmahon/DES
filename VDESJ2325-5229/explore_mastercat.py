from __future__ import (print_function, division)
"""

explore the catalogue data for lensed quasar master catalogue

TODO:
weighted mean of the DES position per waveband

join the mastercat to each catalogue file


"""

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

from librgm.plot_radec import plot_radec
from librgm.xmatch import xmatch_cat
from librgm.xmatch import xmatch_checkplot
from librgm.xmatch import xmatch_checkplots
from librgm.xmatch import xmatch_selfcheck


from plot_radec_descat import *
from plot_radec_vhscat import *
from plot_radec_wisecat import *
from plot_radec_gaiacat import *


def getargs():
    """
    parse command line arguements

    """

    import argparse

    description = ''
    epilog = ''
    parser =  argparse.ArgumentParser(
        description=description, epilog=epilog,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    configfile_default = 'explore_mastercat.cfg'
    parser.add_argument("--configfile",
                        default=configfile_default,
                        help="configuration file")

    format_default = 'COADD'
    parser.add_argument("--format", default=format_default,
        help="format as SE or COADD [Default]")

    release_default = 'Y1A1'
    parser.add_argument("--release", default=release_default,
        help="release  as a string e.g. Y1A1 [Default]")

    parser.set_defaults(source='VDESJ2325-5229')
    parser.add_argument("--source",
                        help="source  as a string")

    parser.set_defaults(summaryfile=None)
    parser.add_argument("--summaryfile",
                        help="summary file with sample")

    parser.set_defaults(mastercat=None)
    parser.add_argument("--mastercat",
                        help="Analysis of mastercat")

    parser.set_defaults(band='i')
    parser.add_argument("--band",
                        help="waveband [g, r, i, z, Y]")

    parser.set_defaults(wise=False)
    parser.add_argument("--wise", action='store_true',
        help="wise overplot option")

    parser.add_argument("--invert_xaxis", action='store_true',
                        default=False,
                        help="invert x-axis")

    parser.add_argument("--invert_yaxis", action='store_true',
                        default=False,
                        help="invert y-axis")

    parser.set_defaults(cutout=False)
    parser.add_argument("--cutout", action='store_true',
        help="cutout option")

    parser.set_defaults(list=False)
    parser.add_argument("--list", action='store_true',
        help="list the lenses")

    parser.add_argument("--plotimages", action='store_true',
        default=False, help="plot images option")

    parser.add_argument("--showplots",
                        action='store_true',
                        default=False,
                        help="showplots option")

    # 38 pixels is 10 arcsec in DES Y1A1 image with
    # single epoch pixel size 0.267 arcsec/pixel
    parser.add_argument("--size",
                        type=int,
                        default=8,
                        help="Size in arcsecs")

    parser.add_argument("--zoom_size",
                        type=int,
                        default=4,
                        help="Zoom in size in arcsecs")

    parser.add_argument("--size_pixels",
                        type=int,
                        default=38,
                        help="Size in pixels")

    parser.add_argument("--edgedetection",
                        action='store_true',
                        default=False,
                        help="Segmentation map edge detection  option")


    parser.set_defaults(zoom=False)
    parser.add_argument("--zoom", action='store_true',
                        default=False, help="zoom option")

    parser.set_defaults(xkcd=False)
    parser.add_argument("--xkcd", action='store_true',
        dest='xkcd', help="xkcd cartoon plot style")

    parser.set_defaults(debug=False)
    parser.add_argument("--debug",
                        action='store_true',
                        help="debug option")

    parser.set_defaults(verbose=False)
    parser.add_argument("--verbose", action='store_true',
        dest='verbose', help="verbose option")

    args = parser.parse_args()

    return args


def getconfig(configfile=None, debug=False):
    """
    read config file

    """
    import configparser

    from configparser import SafeConfigParser
    # parser = SafeConfigParser()

    if debug:
        print('configfile:', configfile)

    # read the configuration file
    config = configparser.RawConfigParser()
    config.read(configfile)

    if debug:
        print('configfile:', configfile)
        print('sections:', config.sections())
        for section_name in config.sections():
            print('Section:', section_name)
            print('Options:', config.options(section_name))
            for name, value in config.items(section_name):
                print('  %s = %s' % (name, value))
        print()


    return config


if __name__ == "__main__":

    import configparser
    import pprint

    t0 = time.time()

    # TODO
    # plot_compass_arrow(direction='NS', length=2.0)
    # key=raw_input("Enter any key to continue: ")

    args = getargs()

    debug = args.debug
    verbose = args.verbose
    if debug or verbose:
        pprint.pprint(args)
        if debug:
            key = raw_input("Enter any key to continue: ")

    configfile = args.configfile
    config = getconfig(configfile=configfile, debug=debug)

    cutout = args.cutout
    band = args.band
    format = args.format
    plotimages = args.plotimages
    release = args.release
    showplots = args.showplots

    # cutout size in arc seconds
    size = args.size
    # size_zoom = args.zoom_size
    # size = args.size_pixels

    source = args.source
    zoom = args.zoom

    if args.xkcd: plt.xkcd()

    if debug:
        key = raw_input("Enter any key to continue: ")

    xrange = [-size/2.0, size/2.0]
    yrange = [-size/2.0, size/2.0]

    # colours to use for grizY plots; could move to config file
    colors = ['blue','green','orange','red', 'maroon']

    path_mastercat = config.get('MASTERCAT', 'path_mastercat')
    filename_mastercat = config.get('MASTERCAT', 'filename_mastercat')
    infile_mastercat = path_mastercat + '/' + filename_mastercat
    print('Reading mastercat file:', infile_mastercat)
    mastercat = Table.read(infile_mastercat)
    mastercat.info('stats')

    path_descat = config.get('MASTERCAT', 'path_descat')
    filename_descat = config.get('MASTERCAT', 'filename_descat')
    infile_descat = path_descat + '/' + filename_descat
    print('Reading descat file:', infile_descat)
    descat = Table.read(infile_descat)
    descat.info('stats')

    wisedata = descat
    vhsdata = descat

    path_gaiacat = config.get('MASTERCAT', 'path_gaiacat')
    filename_gaiacat = config.get('MASTERCAT', 'filename_gaiacat')
    infile_gaiacat = path_gaiacat + '/' + filename_gaiacat
    print('Reading gaiacat file:', infile_gaiacat)
    gaiacat = Table.read(infile_gaiacat)
    gaiacat.info('stats')

    table1 = descat
    table2 = mastercat

    colnames_radec1 = ['RA', 'DEC']
    colnames_radec2 = ['RA', 'Dec']

    # I do not use the xmatch at the moment since I just trawl through the
    # mastercat sources and descat sources outside the regions get ignored
    # at the plot stage
    if debug:
        help(xmatch_cat)
    print("Elapsed time %.3f seconds" % (time.time() - t0))
    """RA, Dec nearest xmatch for two lists; returns pointers """
    idxmatch, dr = xmatch_cat(table1=table1,
                              table2=table2,
                              colnames_radec1=colnames_radec1,
                              colnames_radec2=colnames_radec2,
                              selfmatch=False,
                              stats=True,
                              debug=debug,
                              verbose=True)
    print("Elapsed time %.3f seconds" % (time.time() - t0))

    dr_median = np.median(dr)
    dr_mad_std = mad_std(dr)
    numpoints = len(dr)
    print(len(dr), dr_median, dr_mad_std)

    itest = np.unique(idxmatch)
    print('Unique idxmatch:', len(itest), len(idxmatch))

    itest = np.unique(table1['row_id'])
    print('Unique row_id:', len(itest))


    for icount, id in enumerate(itest):
        print(icount + 1, id,
              mastercat['RA'][id], mastercat['Dec'][id])
        precision = 1
        print(icount + 1, id,
              Angle(mastercat['RA'][id],
                    unit=u.deg).to_string(unit=u.hour,
                    precision=precision+1, sep=' ', pad=True),
              Angle(mastercat['Dec'][id],
                    unit=u.deg).to_string(unit=u.degree,
                    precision=precision, sep=' ',
                    pad=True, alwayssign=True))

        imatch = (idxmatch == id)
        print(icount + 1, id, len(imatch), len(descat[imatch]))
        print(np.average(descat['RA'][imatch]))
        print(np.average(descat['DEC'][imatch]))

    if debug:
        key = raw_input("Enter any key to continue: ")


    showplot = True
    data = descat
    gaiadata = gaiacat

    for icount, id in enumerate(itest):
        print(icount + 1, id,
              mastercat['RA'][id], mastercat['Dec'][id])

        sourceName = Angle(mastercat['RA'][id],
                           unit=u.deg).to_string(unit=u.hour,
                           precision=2, sep='', pad=True) + \
                     Angle(mastercat['Dec'][id],
                           unit=u.deg).to_string(unit=u.degree,
                           precision=1, sep='',
                           pad=True, alwayssign=True)

        imatch = (idxmatch == id)
        print(icount + 1, id, len(imatch), len(descat[imatch]))

        RA_average = np.average(descat['RA'][imatch])
        DEC_average = np.average(descat['DEC'][imatch])

        radec_centre = (mastercat['RA'][id], mastercat['Dec'][id])
        radec_centre = (RA_average, DEC_average)

        # Note we keep plt handle for overlays; this could cause problems
        plt = plot_radec_descat(data=data,
                                radius=0.45,
                                source=sourceName,
                                radec_centre=radec_centre,
                                xrange=xrange,
                                yrange=yrange,
                                coadd=True,
                                multiBand=True,
                                singleEpoch=False,
                                showplot=False,
                                debug=debug)

        if args.invert_xaxis:
            plt.gca().invert_xaxis()

        if wisedata is not None:
            colnames_radec = ['RA_WISE', 'DEC_WISE']
            plt = plot_radec_wisecat(data=wisedata, radius=3.0,
                                     source=sourceName,
                                     radec_centre=radec_centre,
                                     colnames_radec=colnames_radec,
                                     xrange=xrange,
                                     yrange=yrange,
                                     showplot=False,
                                     plt=plt,
                                     debug=debug)

        if vhsdata is not None:
            colnames_radec = ['RA_VHS', 'DEC_VHS']
            plt = plot_radec_vhscat(data=vhsdata, radius=0.45,
                                    source=sourceName,
                                    radec_centre=radec_centre,
                                    colnames_radec=colnames_radec,
                                    xrange=xrange,
                                    yrange=yrange,
                                    showplot=False,
                                    plt=plt,
                                    debug=debug)

        if gaiadata is not None:
            colnames_radec = ['gaia_ra', 'gaia_dec']
            plt = plot_radec_gaiacat(data=gaiadata, radius=0.10,
                                    source=sourceName,
                                    radec_centre=radec_centre,
                                    colnames_radec=colnames_radec,
                                    xrange=xrange,
                                    yrange=yrange,
                                    showplot=False,
                                    plt=plt,
                                    debug=debug)


        if args.invert_xaxis:
            plt.gca().invert_xaxis()

        plotfile = source + '_COADD_radec.png'
        plt.savefig(plotfile)
        print('Saving: ', plotfile)


        plt.show()
