"""

Compare DES catalogue parameters to investigate which parameters are
different in each waveband

Based on code by Fernanda Ostrovski


"""

from __future__ import print_function, division

import sys

import numpy as np
import pylab as plt
from plotid import plotid

from astropy.table import Table, hstack


def multihist(table=None, colname_prefix=None, tile=None, release=None,
        figname=None, lim=False, thetas=False, frad=False):
    """


    """

    plt.subplots_adjust(bottom=0.2, wspace=0.25)
    plt.title(tile + ' ' + release)

    xmin = np.nanmin(table[colname_prefix + '_G'])
    xmax = np.nanmax(table[colname_prefix + '_G'])
    print('column name prefix:', colname_prefix)
    print('minimum (G):', xmin)
    print('maximum (G):', xmax)
    print('range(G):', xmax - xmin)
    print('mean(G):    ', np.nanmean(table[colname_prefix + '_G']))
    print('median(G):  ', np.nanmedian(table[colname_prefix + '_G']))

    binfactor = 20
    if thetas:
        bins = np.arange(np.nanmin(table[colname_prefix + '_G']),
                         np.nanmax(table[colname_prefix + '_G']),
                         np.abs(np.nanmean(table[colname_prefix + '_G'])))
    else:
        bins = np.arange(np.nanmin(table[colname_prefix + '_G']),
                         np.nanmax(table[colname_prefix + '_G']),
                         np.abs(np.nanmean(table[colname_prefix + '_G'])/binfactor))

    bins = 100
    xdata = table[colname_prefix + '_G']
    xdata = xdata[~np.isnan(xdata)]
    plt.hist(xdata, bins,
             histtype='step', color='#2e64fe', label='g')

    # lazy fix
    t = table
    x = colname_prefix
    # bins=np.arange(min(t[x+'_R']),max(t[x+'_R']),np.mean(t[x+'_G'])/10)
    print('Overplotting r')
    xdata = table[colname_prefix + '_R']
    xdata = xdata[~np.isnan(xdata)]
    plt.hist(xdata, bins,
             histtype='step', color='#92cd00',label='r')
    # bins=np.arange(min(t[x+'_I']),max(t[x+'_I']),np.mean(t[x+'_G'])/10)
    print('Overplotting i')
    xdata = table[colname_prefix + '_I']
    xdata = xdata[~np.isnan(xdata)]
    plt.hist(xdata, bins,
             histtype='step', color='#ffbf00', label='i')
    # bins=np.arange(min(t[x+'_Z']),max(t[x+'_Z']),np.mean(t[x+'_G'])/10)
    print('Overplotting z')
    xdata = table[colname_prefix + '_Z']
    xdata = xdata[~np.isnan(xdata)]
    plt.hist(xdata, bins,
             histtype='step', color='#eb7a3c', label='z')
    # bins=np.arange(min(t[x+'_Y']),max(t[x+'_Y']),np.mean(t[x+'_Y'])/10)
    print('Overplotting Y')
    xdata = table[colname_prefix + '_Y']
    xdata = xdata[~np.isnan(xdata)]
    plt.hist(xdata, bins,
             histtype='step', color='#cc0084',label='Y')
    plt.xlabel(x)
    plotid()

    if lim:
                plt.xlim(min(t[x+'_G']),np.median(t[x+'_G'])+5*(np.std(t[x+'_G'])))
    if frad:
                #plt.xlim(np.median(t[x+'_G'])-(np.std(t[x+'_G'])),np.median(t[x+'_G'])+(np.std(t[x+'_G'])))
                plt.xlim(0.0,10)

    plt.legend(prop={'size':10}, loc='upper center',
               bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=False, ncol=5)

    # plt.show()
    print('Save:', figname)
    plt.savefig(figname)
    plt.clf()


def fxy(tS,tY,x,y,xlab,ylab,figname,yinvert=False, lim=False,frad=False):
        plt.scatter(tS[x],tS[y],s=1,edgecolor='none',c='#2e64fe', label='SVA1')
        plt.scatter(tY[x],tY[y],s=1,edgecolor='none',c='#ff00a5',label='Y1A1')
        plt.xlabel(xlab,fontsize='17')
        plt.ylabel(ylab,fontsize='17')
        px=np.arange(min(tS[x]),max(tS[x])+1,0.1)
        plt.plot(px,px,c='k')
        if lim:
                plt.xlim(min(tS[x]),np.median(tS[x])+5*(np.std(tS[x])))
		plt.ylim(min(tS[x]),np.median(tS[x])+5*(np.std(tS[x])))
        if frad:
                #plt.xlim(np.median(t[x+'_G'])-(np.std(t[x+'_G'])),np.median(t[x+'_G'])+(np.std(t[x+'_G'])))
                plt.xlim(1.0,12)
		plt.ylim(1.0,12)
	if yinvert:
		plt.ylim(14,26)
                plt.gca().invert_yaxis()
        plt.title(tile)
        plotid()
        plt.legend(prop={'size':10})
        plt.show()
        #plt.savefig(figname)
        plt.clf()

def abtheta(x2,y2,xy):
        a2=((x2+y2)/2)+np.sqrt(((x2-y2)/2)*((x2-y2)/2)+xy*xy)
        b2=((x2+y2)/2)-np.sqrt(((x2-y2)/2)*((x2-y2)/2)+xy*xy)
        a=np.sqrt(a2)
        b=np.sqrt(b2)
        ellip=1-(b/a)
        return a, b, ellip


def des_joinbands(tile=None, release=None, suffix='r2587p01'):

        infile = '/data/desardata/'+release+'/'+tile+'/'+tile+'_g_cat.fits'
        if release == 'Y3A1':
            infile = '/data/desardata/' + release + '/' + tile + '/' + tile+ '_' + suffix + '_g_cat.fits'
        g = Table.read(infile)

        infile = '/data/desardata/'+release+'/'+tile+'/'+tile+'_r_cat.fits'
        if release == 'Y3A1':
            infile = '/data/desardata/' + release + '/' + tile + '/' + tile+ '_' + suffix + '_r_cat.fits'
        r = Table.read(infile)

        infile = '/data/desardata/'+release+'/'+tile+'/'+tile+'_i_cat.fits'
        if release == 'Y3A1':
            infile = '/data/desardata/' + release + '/' + tile + '/' + tile+ '_' + suffix + '_i_cat.fits'
        i = Table.read(infile)

        infile = '/data/desardata/'+release+'/'+tile+'/'+tile+'_z_cat.fits'
        if release == 'Y3A1':
            infile = '/data/desardata/' + release + '/' + tile + '/' + tile+ '_' + suffix + '_z_cat.fits'
        z = Table.read(infile)

        infile = '/data/desardata/'+release+'/'+tile+'/'+tile+'_Y_cat.fits'
        if release == 'Y3A1':
            infile = '/data/desardata/' + release + '/' + tile + '/' + tile+ '_' + suffix + '_Y_cat.fits'
        Y = Table.read(infile)

        # stack all the columns across the DES wavebands
        t = hstack([g,r,i,z,Y], table_names=['G','R','I','Z','Y'])

        t.info()

        t.write('/data/desardata/'+release+'/'+tile+'/'+tile+'_merged_cat.fits', overwrite=True)

        print(release, tile, 'done')


def get_colname_prefixes(table=None):
    """
    return the unique list of colname prefixes

    """
    colnames = table.colnames
    print(type(colnames), len(colnames))
    print(colnames[0], colnames[-1])
    prefixes = [word[:-2] for word in colnames]
    print(prefixes[0], prefixes[-1])

    print(len(prefixes))
    prefixes = np.asarray(prefixes)
    print(len(prefixes))

    unique_prefixes = np.unique(prefixes)
    print(unique_prefixes[0], unique_prefixes[-1])
    print(len(unique_prefixes))

    return unique_prefixes

if __name__ == '__main__':
    """
    compare SVA1 and Y1A1 catalogue for the same tile

    """

    # tile='DES0453-4457'
    # tile='DES0449-4706'
    tile = 'DES2327-5248'
    dirY = '/data/desardata/Y1A1/' + tile + '/'
    dirY = '/data/desardata/Y3A1/' + tile + '/'
    dirS = '/data/desardata/SVA1/' + tile + '/'

    releases=['SVA1','Y1A1']
    releases=['Y1A1']
    releases=['SVA1']
    releases=['Y3A1']
    for release in releases:
        des_joinbands(tile=tile, release=release, suffix='r2587p01')

    DEBUG = True
    infile =  dirY + tile + '_merged_cat.fits'
    print('Read:', infile)
    tY = Table.read(infile)
    colname_prefixes = get_colname_prefixes(table=tY)
    if DEBUG:
        raw_input("Enter any key to continue: ")

    tS = []
    # tS = Table.read(dirS + tile+'_merged_cat.fits')

    tY = tY[tY['MAG_AUTO_I']<21.0]
    # tS = tS[tS['MAG_AUTO_I']<19.0]

    tS = []
    print(len(tS), len(tY))


    # histogram of flux_radius
    # multhist(tS,'FLUX_RADIUS','SVA1',tile+'_fluxradius_SVA1.png',lim=False,thetas=False,frad=True)
    release = 'Y3A1'
    for colname_prefix in colname_prefixes:
        figname = tile + '_' + colname_prefix + '_' + release + '_multihist' + '.png'
        multihist(table=tY, colname_prefix=colname_prefix, release=release, tile=tile,
                  figname=figname,
                  lim=False, thetas=False, frad=False)

    sys.exit()

    # mag x flux_radius (ratios?)
    # fxy(tS,tY,'FLUX_RADIUS_G','MAG_AUTO_G','FLUX_RADIUS_G','MAG_AUTO_G',tile+'_fluxradius_mag_g.png',yinvert=True, lim=False,frad=True)
    # fxy(tS,tY,'FLUX_RADIUS_I','MAG_AUTO_I','FLUX_RADIUS_I','MAG_AUTO_I',tile+'_fluxradius_mag_i.png',yinvert=True, lim=False,frad=True)

    #a/b/fluxradius g x i
    bands=['G','R','I','Z','Y']
    for band in bands:
        x2=np.array(tY['X2WIN_IMAGE_' + band])
        y2=np.array(tY['Y2WIN_IMAGE_' + band])
        xy=np.array(tY['XYWIN_IMAGE_' + band])
        a,b,ellip=abtheta(x2,y2,xy)
        tY['AWIN_IMAGE_'+band]=np.array(a)
        tY['BWIN_IMAGE_'+band]=np.array(b)
        tY['ELLIPTICTYWIN_'+band]=np.array(ellip)

    # fxy(tS,tY,'AWIN_IMAGE_G','AWIN_IMAGE_I','AWIN_IMAGE_G','AWIN_IMAGE_I',tile+'_awin_gi.png',yinvert=False, lim=True,frad=False)
    # fxy(tS,tY,'BWIN_IMAGE_G','BWIN_IMAGE_I','BWIN_IMAGE_G','BWIN_IMAGE_I',tile+'_bwin_gi.png',yinvert=False, lim=True,frad=False)
    fxy(tS,tY,'FLUX_RADIUS_G','FLUX_RADIUS_I','FLUX_RADIUS_G','FLUX_RADIUS_I',tile+'_fluxradius_gi.png',yinvert=False, lim=False,frad=True)

    # mag x a,b
    # fxy(tS,tY,'AWIN_IMAGE_G','MAG_AUTO_G','AWIN_IMAGE_G','MAG_AUTO_G',tile+'_awin_mag_g.png',yinvert=True, lim=True,frad=True)
    # fxy(tS,tY,'AWIN_IMAGE_I','MAG_AUTO_I','AWIN_IMAGE_I','MAG_AUTO_I',tile+'_awin_mag_i.png',yinvert=True, lim=True,frad=True)
    # fxy(tS,tY,'BWIN_IMAGE_G','MAG_AUTO_G','BWIN_IMAGE_G','MAG_AUTO_G',tile+'_bwin_mag_g.png',yinvert=True, lim=True,frad=True)
    # fxy(tS,tY,'BWIN_IMAGE_I','MAG_AUTO_I','BWIN_IMAGE_I','MAG_AUTO_I',tile+'_bwin_mag_i.png',yinvert=True, lim=True,frad=True)
