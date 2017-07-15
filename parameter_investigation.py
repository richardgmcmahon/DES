from __future__ import (print_function, division)

import numpy as np
import matplotlib.pyplot as plt

from astropy.table import Table
import astropy.io.fits as fits
import astropy.io.fits.compression

import Possibles_Analysis as PA
import plotid
import stats
from match_lists import match_lists as ml


def make_hist(xs, col, units, comment, band, file_start, out_path,
              zoom = False, save = True):

    fig = plt.figure()
    ids = np.where( (xs == xs) )[0]
    xs = xs[ids]
    pers = np.percentile(xs, [1.0, 99.0])
    keeps = np.where( (xs < pers[1]) & (xs > pers[0]))[0]
    if zoom and len(keeps) > 1:
        xs1 = xs[keeps]
        nper = len(xs1)
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122, sharey = ax1)
        ax2.get_yaxis().set_visible(False)
        ax2.hist(xs1, bins = 100, log = True, range = (min(xs1), max(xs1)))
        ax2.set_title("1st - 99th %tile: " + str(nper))
        labels2 = ax2.get_xticks()
        ax2.set_xticklabels(labels2, rotation = 270)
    else:
        ax1 = fig.add_subplot(111)

    nr = len(xs)
    ax1.hist(xs, bins = 100, log = True, range = (min(xs), max(xs)))
    labels1 = ax1.get_xticks()[:-1]
    ax1.set_xticklabels(labels1, rotation = 270)
    ax1.text(0.2, 0.7, "Min: " + str(min(xs)) + "\nMax: " + str(max(xs)) +  "\nMedian: " + str(np.median(xs)) + "\nSigma MAD: " + str(1.4826*stats.MAD(xs, np.median(xs))) + "\n1st %ile: " + str(pers[0]) + "\n99th %ile: " + str(pers[1]), transform=ax1.transAxes, bbox=dict(facecolor='blue', alpha=0.2))
    ax1.set_title("All points: " + str(nr))
    ax1.text(0.5, 0.05, col + " / " + units + "\n" + comment, ha = "center", transform = fig.transFigure)
    ax1.set_ylabel("Frequency")
    print col, file_start, band
    fig.suptitle(col + " " + file_start + band)
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2, wspace = 0.0)
    plotid.plotid()
    fig = plt.gcf()
    fig.set_size_inches(10.0,8.0)
    if save:
        plt.savefig(out_path + col + "_" + file_start + band + ".png")
        plt.close()
    else:
        plt.show()

def histograms(dir, file_start, file_end, cols, bands, out_path,
               zoom = False, save = True):

    for band in bands:
        t = Table.read(dir + file_start + band + file_end)
        hdr = fits.open(dir + file_start + band + file_end)

        for col in cols:
            column_list = list(t.columns)
            print col, band

            if "-" not in col:
                xs = np.array(t[col], dtype = np.float64)
                units = str(t[col].units)
                i = column_list.index(col) + 1
                comment = "(" + hdr[1].header.comments["TTYPE" + str(i)] + ")"
                if len(xs.shape) > 1:
                    n = 0
                    while n < len(t[col][0]):
                        xs = t[col][:,n]
                        col1 = col + "_" + str(n+1)
                        n += 1
                        print col1
                        make_hist(xs, col1, units, comment, band, file_start + file_end, out_path, zoom = zoom, save = save)
                else:
                    make_hist(xs, col, units, comment, band, file_start + file_end, out_path, zoom = zoom, save = save)

            else:
                l = col.index("-")
                col1 = col[:l]
                col2 = col[l+1:]
                xs = t[col1] - t[col2]
                units = str(t[col1].units)
                i1 = column_list.index(col1) + 1
                i2 = column_list.index(col2) + 1
                comment = "(" + hdr[1].header.comments["TTYPE" + str(i1)] + " " + hdr[1].header.comments["TTYPE" + str(i2)] + ")"
                if len(xs) > 1:
                    #for xs in t[col]:
                    make_hist(xs, col, units, comment, band, file_start + file_end, out_path, zoom = zoom, save = save)
                else:
                    print "No data"
                    #make_hist(xs, col, units, comment, band, file_start + file_end, out_path, zoom = zoom, save = save)

def AB_image(dir, file_start, file_end, bands):

    for band in bands:
        with fits.open(dir + file_start + band + file_end) as hlist:
            data = hlist[1].data
            if band == "g":
                A0s = data["A_IMAGE"]
                B0s = data["B_IMAGE"]
                fwhm0s = data["FWHM_IMAGE"]

            else:
                A1s = data["A_IMAGE"]
                B1s = data["B_IMAGE"]
                fwhm1s = data["FWHM_IMAGE"]
                #print sum(A0s - A1s)
                #print sum(B0s - B1s)
                #print fwhm1s**2/(A1s**2 + B1s**2)
                #plt.plot(fwhm1s**2, (A1s**2 + B1s**2), "k.")
                #plt.plot(fwhm1s**2, ((A1s+B1s)/2.0)**2, "r.")
                isos = data["ISOAREA_IMAGE"]
                plt.plot(A1s*B1s, isos, "k.", ms = 1)
                plt.xlabel("A_IMAGE * B_IMAGE")
                plt.ylabel("ISOAREA_IMAGE")
                plt.title("Made with parameter_investigation.py for the " + band + " band")
                plt.show()

            #if band == "i":
                #for (n,A) in enumerate(A1s):
                #    if A > 6.0:
                #        fig = PA.cutout_image("", data["ALPHAWIN_J2000"][n], data["DELTAWIN_J2000"][n], "DES0332-2749", "20130305000001_DES0332-2749", data["NUMBER"][n], cat_info = True)
                #        plt.show()
                #for (n,B) in enumerate(B1s):
                #    if B > 4.0:
                #        fig = PA.cutout_image("", data["ALPHAWIN_J2000"][n], data["DELTAWIN_J2000"][n], "DES0332-2749", "20130305000001_DES0332-2749",     data["NUMBER"][n], cat_info = True)
                #        plt.show()

def kron_radius(dir, file_start, file_end, bands):

    for band in bands:
        with fits.open(dir + file_start + band + file_end) as hlist:
            data = hlist[1].data
            ks = data["KRON_RADIUS"]
            print min(ks), max(ks)
            if band == "g":
                k0s = ks
            else:
                k1s = ks
                diffs = (k0s - k1s)
            ids = np.where( (ks == 0) )[0]
            ras = data["ALPHAWIN_J2000"][ids]
            decs = data["DELTAWIN_J2000"][ids]
            tile = "DES0332-2749"
            run = "20130305000001_DES0332-2749"
            plt.plot(ras, decs, "k.")
            plt.show()

            n = 1000
            while n < 1010:
                RA = ras[n]
                DEC = decs[n]
                id = data["NUMBER"][n]
                PA.cutout_image("", RA, DEC, tile, run, id, save = False, cat_info = True)
                plt.show()
                n += 1

def elongation(dir, file_start, file_end, bands):

    band = "i"
    with fits.open(dir + file_start + band + file_end) as hlist:
        data = hlist[1].data
        es = data["ELONGATION"]

    for (n,e) in enumerate(es):
        if e > 3.0:
            fig = PA.cutout_image("", data["ALPHAWIN_J2000"][n], data["DELTAWIN_J2000"][n], "DES0332-2749", "20130305000001_DES0332-2749",     data["NUMBER"][n], cat_info = True)
            plt.show()

def XY_min_max(dir, file_start, file_end, bands):
    band = "i"
    with fits.open(dir + file_start + band + file_end) as hlist:
        data = hlist[1].data
        xs = data["XMAX_IMAGE"]-data["XMIN_IMAGE"]

    for (n,x) in enumerate(xs):
        if x > 50.0:
            fig = PA.cutout_image("", data["ALPHAWIN_J2000"][n], data["DELTAWIN_J2000"][n], "DES0332-2749", "20130305000001_DES0332-2749",     data["NUMBER"][n], cat_info = True)
            plt.show()

def isoarea(dir, file_start, file_end, bands):
    band = "i"
    with fits.open(dir + file_start + band + file_end) as hlist:
        data = hlist[1].data
        xs = data["ISOAREA_IMAGE"]

    for (n,x) in enumerate(xs):
        if x == 0:
            fig = PA.cutout_image("", data["ALPHAWIN_J2000"][n], data["DELTAWIN_J2000"][n], "DES0332-2749", "20130305000001_DES0332-2749", data["NUMBER"][n], cat_info = True)
            #fig = PA.cutout_zoom("", data["ALPHAWIN_J2000"][n], data["DELTAWIN_J2000"][n], "DES0332-2749", data["NUMBER    "][n], save = False)
            plt.show()

def petro_radius(dir, file_start, file_end, bands):
    band = "i"
    with fits.open(dir + file_start + band + file_end) as hlist:
        data = hlist[1].data
        xs = data["PETRO_RADIUS"]

    for (n,x) in enumerate(xs):
        if x > 10:
            fig = PA.cutout_image("", data["ALPHAWIN_J2000"][n], data["DELTAWIN_J2000"][n], "DES0332-2749", "20130305000001_DES0332-2749", data["NUMBER"][n], cat_info = True)
            #fig = PA.cutout_zoom("", data["ALPHAWIN_J2000"][n], data["DELTAWIN_J2000"][n], "DES0332-2749", data["NUMBER    "][n], save = False)
            plt.show()

def FWHM(dir, file_start, file_end, bands):
    band = "i"
    with fits.open(dir + file_start + band + file_end) as hlist:
        data = hlist[1].data
        xs = data["FWHM_IMAGE"]

    for (n,x) in enumerate(xs):
        if x == 0:
            fig = PA.cutout_image("", data["ALPHAWIN_J2000"][n], data["DELTAWIN_J2000"][n], "DES0332-2749", "20130305000001_DES0332-2749", data["NUMBER"][n], cat_info = True)
            #fig = PA.cutout_zoom("", data["ALPHAWIN_J2000"][n], data["DELTAWIN_J2000"][n], "DES0332-2749", data["NUMBER    "][n], save = False)
            plt.show()

def flux_radius(dir, file_start, file_end, bands):
    band = "i"
    with fits.open(dir + file_start + band + file_end) as hlist:
        data = hlist[1].data
        xs = data["FLUX_RADIUS"]

    for (n,x) in enumerate(xs):
        if x < 0:
            fig = PA.cutout_image("", data["ALPHAWIN_J2000"][n], data["DELTAWIN_J2000"][n], "DES0332-2749", "20130305000001_DES0332-2749", data["NUMBER"][n], cat_info = True)
            #fig = PA.cutout_zoom("", data["ALPHAWIN_J2000"][n], data["DELTAWIN_J2000"][n], "DES0332-2749", data["NUMBER    "][n], save = False)
            plt.show()

def fluxes(dir, file_start, file_end, bands):
    for band in bands:
        with fits.open(dir + file_start + band + file_end) as hlist:
            data = hlist[1].data
            fauto = data["FLUX_AUTO"]
            fmax = data["FLUX_MAX"]
            fmodel = data["FLUX_MODEL"]
            fpsf = data["FLUX_PSF"]
            mmax = data["MU_MAX"]
            thres = data["THRESHOLD"]
            mthres = data["MU_THRESHOLD"]

            #xs = np.log10(fmodel)
            xs = np.log10(fpsf)
            #xs = np.log10(fauto)
            #ys = np.log10(fmodel/fmax)
            ys = np.log10(fmodel/fmax)
            #ys = np.log10(fauto/fmax)
            #xlabel = "Log10(FLUX_MODEL)"
            xlabel = "Log10(FLUX_PSF)"
            #xlabel = "Log10(FLUX_AUTO)"
            #ylabel = "log10(FLUX_MODEL_div_FLUX_MAX)"
            ylabel = "log10(FLUX_PSF_div_FLUX_MAX)"
            #ylabel = "Log10(FLUX_AUTO_div_FLUX_MAX)"

            """
            fig = plt.figure()
            ax1 = fig.add_subplot(311)
            ax1.plot(xs, ys, "k.", ms = 1)
            ax1.axes.get_xaxis().set_visible(False)
            ax2 = fig.add_subplot(312, sharex = ax1)
            ax2.plot(xs, ys, "k.", ms = 1)
            labels2 = ax2.get_yticks()[:-1]
            ax2.set_yticklabels(labels2)
            pers = np.percentile(ys, [0.5, 99.5])
            ax2.set_ylim(pers[0], pers[1])
            ax2.axes.get_xaxis().set_visible(False)
            ax3 = fig.add_subplot(313, sharex = ax1)
            ax3.plot(xs, ys, "k.", ms = 1)
            pers = np.percentile(ys, [2.0, 98.0])
            ax3.set_ylim(pers[0], pers[1])
            ax3.set_xlabel(xlabel)
            labels3 = ax3.get_yticks()[:-1]
            ax3.set_yticklabels(labels3)
            plotid.plotid()
            plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.1, wspace = 0.0, hspace = 0.0)
            fig = plt.gcf()
            fig.set_size_inches(10.0,20.0)
            ax2.text(0.05, 0.5, ylabel, va = "center", transform = fig.transFigure, rotation = "vertical")
            plt.suptitle(file_start + band)
            #plt.savefig("/home/sr525/Graphs/Parameters/" + xlabel + "_v_" + ylabel + "_" + file_start + band + ".png")
            #plt.close()
            plt.show()
            """

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(xs, ys, "k.", ms = 1)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title(file_start + band)
            plotid.plotid()
            medx = np.median(xs)
            medy = np.median(ys)
            MADx = stats.MAD(xs, medx)
            MADy = stats.MAD(ys, medy)
            maxx = max(xs)
            minx = min(xs)
            maxy = max(ys)
            miny = min(ys)
            ax.text(0.1, 0.7, "Range x: %0.2f " % (min(xs)) + "to %0.2f\n" % (max(xs)) + "Range y: %0.2f " % (min(ys)) + "- %0.2f\n" % (max(ys)) + "Medians x: %0.2f" % (medx) + " y: %0.2f\n" % (medy) + "MADs x: %0.2f" % (MADx) + " y: %0.2f\n" % (MADy) + "Sigma MADs x: %0.2f" % (1.4826*MADx) + " y: %0.2f" % (1.4826*MADy), transform=ax.transAxes, bbox=dict(facecolor='black', alpha=0.2))
            plt.savefig("/home/sr525/Graphs/Parameters/" + xlabel + "_v_" + ylabel + "_" + file_start + band + ".png")
            plt.close()
            #plt.show()

            """
            plt.plot(np.log10(fmax), fmax/mmax, "k.", ms = 1)
            plt.xlabel("log10(FLUX_MAX)")
            plt.ylabel("FLUX_MAX / MU_MAX")
            plt.title(file_start + band)
            plt.savefig("/home/sr525/Graphs/Parameters/Log10(FLUX_MAX)_v_FLUX_MAX_div_MU_MAX.png")
            plt.close()
            #plt.show()

            plt.plot((thres), thres/mthres, "k.", ms = 1)
            plt.xlabel("THRESHOLD")
            plt.ylabel("THRESHOLD / MU_THRESHOLD")
            plt.title(file_start + band)
            plt.savefig("/home/sr525/Graphs/Parameters/THRESHOLD_v_THRESHOLD_div_MU_THRESHOLD.png")
            plt.close()
            #plt.show()
            """

def ra_dec(dir, file_start, file_end, bands):
    """

    """
    rass = [[],[],[],[],[]]
    decss = [[], [], [], [] ,[]]
    numss = [[], [], [], [], []]
    for band in bands:
        b = bands.index(band)
        with fits.open(dir + file_start + band + file_end) as hlist:
            data = hlist[1].data
            ras = data["ALPHAWIN_J2000"]
            decs = data["DELTAWIN_J2000"]
            print len(ras)
            rass[b] = ras
            decss[b] = decs
            numss[b] = data["NUMBER"]

    t = Table.read(dir + file_start[:-1] + ".fits")
    #dists, inds = ml(rass[0], decss[0], rass[1], decss[1], 0.0006)
    #dists1, inds1 = ml(rass[0], decss[0], rass[2], decss[2], 0.0006)
    #ids = np.where( (inds <> inds1) & (inds < len(rass[0])) )[0]
    #print ids
    #print len(ids)
    #inds = inds[ids]
    dists, inds = ml(t["ALPHAWIN_J2000_G"], t["DELTAWIN_J2000_G"], t["ALPHAWIN_J2000_R"], t["DELTAWIN_J2000_R"], 0.0006)
    ids = np.where( (inds == len(t)) )[0]
    print len(ids)
    t_odd = t[ids]
    n = 0
    #while n < len(t_odd):
    #    fig = PA.cutout_image("", t_odd["ALPHAWIN_J2000_G"][n], t_odd["DELTAWIN_J2000_G"][n], "DES0332-2749", "20130305000001_DES0332-2749", t_odd["COADD_OBJECTS_ID"][n], cat_info = True)
    #    plt.show()
    #    n += 1
    n = 0
    plt.plot(t["RA"], t["DEC"], "k.", ms = 1)
    plt.show()

    #while n < len(t):
        #print t["RA"][n], np.mean([t["ALPHAWIN_J2000_G"][n], t["ALPHAWIN_J2000_R"][n], t["ALPHAWIN_J2000_I"][n], t["ALPHAWIN_J2000_Z"][n]])
        #if t["ALPHAWIN_J2000_Y"][n] > 53.5:
        #    PA.cutout_image("", t["ALPHAWIN_J2000_Z"][n], t["DELTAWIN_J2000_Z"][n], "DES0332-2749", "20130305000001_DES0332-2749", t["COADD_OBJECTS_ID"][n], cat_info = True)
        #    plt.show()
        #n += 1


    r = 0
    for ras in rass:
        print ras[ids][0:10], decss[r][ids][0:10]
        plt.plot(ras, decss[r], "k.", ms = 1)
        plt.show()
        r += 1


def background(dir, file_start, file_end, bands):
    for band in bands:
        with fits.open(dir + file_start + band + ".fits.fz") as hlist:
            im = hlist[1].data
            print np.median(im)

        with fits.open(dir + file_start + band + file_end) as hlist1:
            data = hlist1[1].data
            print np.median(data["BACKGROUND"])


def chi(dir, file_start, file_end, bands):
    for band in bands:
        with fits.open(dir + file_start + band + file_end) as hlist:
            data = hlist[1].data
            xs = data["CHI2_PSF"]

        for (n, x) in enumerate(xs):
            if x > 7e+23:
                fig = PA.cutout_image("", data["ALPHAWIN_J2000"][n],
                                      data["DELTAWIN_J2000"][n],
                                      "DES0332-2749",
                                      "20130305000001_DES0332-2749",
                                      data["NUMBER"][n], cat_info=True)
                #fig = PA.cutout_zoom("", data["ALPHAWIN_J2000"][n],
                #data["DELTAWIN_J2000"][n], "DES0332-2749",
                #data["NUMBER"][n], save = False)

                plt.show()


if __name__ == '__main__':
    """

    """
    #dir = "/data/desardata/SVA1/DES1000+0209/"
    #t = Table.read(dir + "DES1000+0209_i_cat.fits")
    #dir = "/data/desardata/SVA1/DES0332-2749/"
    #cols = t.columns
    #cols = ["A_IMAGE", "B_IMAGE", "XMAX_IMAGE-XMIN_IMAGE", "
    # YMAX_IMAGE-YMIN_IMAGE", "A_WORLD", "B_WORLD", "ISOAREA_WORLD",
    # "ISOAREA_IMAGE", "ISOAREAF_IMAGE", "KRON_RADIUS", "PETRO_RADIUS",
    # "FWHM_IMAGE", "FWHM_WORLD", "THETA_IMAGE", "ELONGATION",
    # "FLUX_RADIUS", "AMODEL_IMAGE", "BMODEL_IMAGE", "THETAMODEL_IMAGE",
    # "AMODEL_WORLD", "BMODEL_WORLD", "THETAMODEL_J2000", "THETA_J2000",
    #"ISOAREAF_WORLD"]
    #cols = ["BMODEL_WORLD"]
    #bands = ["g", "r", "i", "z", "Y"]
    #file_start = "DES0332-2749_"
    #file_end = "_cat.fits"

    #file_start = "DES0332-2749_"
    #file_start = "DES1000+0209_"
    #dir = "/data/desardata/SVA1/" + file_start[:-1] + "/"
    #for file_start in ["DES0332-2749_", "DES1000+0209_", "DES0453-4457_"]:
    #    t = Table.read("/data/desardata/SVA1/" + file_start[:-1] +
    #  "/" + file_start + "i" + file_end)
    #    cols = t.columns
    #    histograms("/data/desardata/SVA1/" + file_start[:-1] +
    # "/", file_start, file_end, cols, bands, zoom = True)
    #AB_image(dir, file_start, file_end, bands)
    #kron_radius(dir, file_start, file_end, bands)
    #elongation(dir, file_start, file_end, bands)
    #XY_min_max(dir, file_start, file_end, bands)
    #isoarea(dir, file_start, file_end, bands)
    #petro_radius(dir, file_start, file_end, bands)
    #FWHM(dir, file_start, file_end, bands)
    #fluxes(dir, file_start, file_end, bands)
    #flux_radius(dir, file_start, file_end, bands)
    #ra_dec(dir, file_start, file_end, bands)
    #background(dir, file_start, file_end, bands)
    #chi(dir, file_start, file_end, bands)
