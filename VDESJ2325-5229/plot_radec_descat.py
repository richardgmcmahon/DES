import sys
import numpy as np

sys.path.append('/home/rgm/soft/python/lib/')
from librgm.plotid import plotid


def moments2ellipse(X2, Y2, XY):
    """
    convert Sextractor 2nd order moments to ellipse parameters

    """

    A2 = (0.5 * (X2 + Y2)) + np.sqrt(np.power((X2 - Y2), 2) + np.power(XY, 2))
    a = np.sqrt(A2)

    B2 = (0.5 * (X2 + Y2)) - np.sqrt(np.power((X2 - Y2), 2) + np.power(XY, 2))
    b = np.sqrt(B2)

    tan2theta = 2.0 * XY / (X2 - Y2)

    # using arctan2 and not arctan
    theta = np.rad2deg(np.arctan2(tan2theta) / 2.0)

    return a, b, theta


def des_parameter_analysis(data=None, release='Y1A1', index=None):
    """
    do some analysis of the shape parameters


    """
    itest = index
        # primary key in table
    COADD_OBJECTS_ID = data['COADD_OBJECTS_ID'][itest]

    # waveband used for the 'master' ra, dec in table
    RADEC_SOURCE_BAND = data['RADEC_SOURCE_BAND'][itest]

    # Detection image shape measurements
    A_IMAGE = data['A_IMAGE'][itest]*PIXEL_SIZE
    B_IMAGE = data['B_IMAGE'][itest]*PIXEL_SIZE
    THETA_IMAGE = data['THETA_IMAGE'][itest]

    KRON_RADIUS = data['KRON_RADIUS'][itest]*PIXEL_SIZE

    XMIN_IMAGE=data['XMIN_IMAGE'][itest]
    XMAX_IMAGE=data['XMAX_IMAGE'][itest]
    YMIN_IMAGE=data['YMIN_IMAGE'][itest]
    YMAX_IMAGE=data['YMAX_IMAGE'][itest]

    XWIN_IMAGE = data['XWIN_IMAGE_' + WAVEBAND][itest]
    YWIN_IMAGE = data['YWIN_IMAGE_' + WAVEBAND][itest]

    #
    NEPOCHS = data['NEPOCHS_' + WAVEBAND][itest]
    MAG_AUTO = data['MAG_AUTO_' + WAVEBAND][itest]
    FLUX_RADIUS = data['FLUX_RADIUS_'+ WAVEBAND][itest] * PIXEL_SIZE

    # shape measurements per waveband
    AWIN_IMAGE = data['AWIN_IMAGE_' + WAVEBAND][itest] * PIXEL_SIZE
    BWIN_IMAGE = data['BWIN_IMAGE_' + WAVEBAND][itest] * PIXEL_SIZE
    THETAWIN_IMAGE = data['THETAWIN_IMAGE_'+ WAVEBAND][itest]

    ELLIP1MODEL_WORLD = data['ELLIP1MODEL_WORLD_' + WAVEBAND][itest]
    ELLIP2MODEL_WORLD = data['ELLIP2MODEL_WORLD_' + WAVEBAND][itest]

    FWHM_WORLD = data['FWHM_WORLD_' + WAVEBAND][itest]

    ISOAREA_WORLD = data['ISOAREA_WORLD_' + WAVEBAND][itest]

    FLAGS = data['FLAGS_' + WAVEBAND][itest]

    return


def plot_radec_descat(data=None, release='Y1A1',
                      sourceName=None,
                      radius = 0.45, alpha=0.2,
                      plot_ellipses=True,
                      line_style=None,
                      colnames_radec=['RA', 'DEC'],
                      radec_centre=None,
                      xrange=[-2.5, 2.5],
                      yrange=[-2.5, 2.5],
                      wavebands=['g', 'r', 'i', 'z', 'y'],
                      colors=['blue', 'green', 'orange', 'red', 'maroon'],
                      coadd=True,
                      multiBand=True,
                      singleEpoch=False,
                      showplot=False,
                      explore_shapepars=False,
                      debug=True):
    """

    plot the DES catalogue sources. The DES catalogues use Sextractor for
    source detection and measurement.

    """

    import numpy as np

    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, Ellipse

    infile = None
    if 'filename' in data.meta:
        infile = data.meta['filename']

    print('Number of rows:', len(data))
    if debug:
        data.info()
        data.info('stats')
        print('data.meta:', data.meta)
        print('filename:', infile)
        print('ra_centre', radec_centre[0])
        print('dec_centre', radec_centre[1])
        print('xrange:', xrange)
        print('yrange:', yrange)

    itest = np.unique(data['TILENAME'])
    print('Tiles:', itest)

    ra = data[colnames_radec[0]]
    dec = data[colnames_radec[1]]

    # this will populate OBJECT_NUMBER and noit crash
    # OBJECT_NUMBER
    # primary key in each tile catalogue
    # OBJECT_NUMBER = data['OBJECT_NUMBER'][itest]
    index_column = -1
    try:
        index_column = data.index_column('OBJECT_NUMBER')
    except:
        pass
    if debug:
        print('index_column:', index_column)
    if index_column > -1:
        OBJECT_NUMBER = data['OBJECT_NUMBER']
    if index_column < 0:
        OBJECT_NUMBER = np.full(len(data), -1)

    COADD_OBJECTS_ID = data['COADD_OBJECTS_ID']

    ra_min = np.min(ra)
    ra_max = np.max(ra)
    # convert to arc seconds and deal with cosine dec
    delta_ra = (ra - radec_centre[0]) * 3600.0 * \
               np.cos(np.deg2rad(radec_centre[1]))

    if debug:
        print('ra_min:', ra_min)
        print('ra_max:', ra_max)
        print('Delta RA range:',
              np.min(delta_ra), np.max(delta_ra))

    dec_min = np.min(dec)
    dec_max = np.max(dec)
    delta_dec = (dec - radec_centre[1]) * 3600.0

    if debug:
        print('dec_min:', dec_min)
        print('dec_max:', dec_max)
        print('Delta Dec range:',
              np.min(delta_dec), np.max(delta_dec))

    if debug:
        key=raw_input("Enter any key to continue: ")

    xdata = delta_ra
    ydata = delta_dec

    print(xrange)
    print(yrange)
    itest = (xdata > xrange[0]) & (xdata < xrange[1]) & \
            (ydata > yrange[0]) & (ydata < yrange[1])

    xdata = xdata[itest]
    ydata = ydata[itest]
    ndata = len(xdata)

    plt.figure(figsize=(8,8))
    plt.axes().set_aspect('equal')
    plt.plot(xdata, ydata, '.k', label='COADD: '+str(ndata))

    plt.xlim(xrange)
    plt.ylim(yrange)


    plt.title(infile, fontsize='medium')
    plt.grid(True)

    WAVEBANDS = ['G','R','I','Z','Y']

    print("id, object_id, WAVEBAND, NEPOCHS, MAG_AUTO, " + \
              "Dra, Ddec, " + \
              "FLUX_RADIUS, KRON_RADIUS, " + \
              "A, B, PA, Aspect Ratio(A/B)")

    # used to convert pixels to arc seconds
    PIXEL_SIZE = 0.27

    # loop through the wavebands
    if debug:
        data.info('stats')
    for iband, WAVEBAND in enumerate(WAVEBANDS):

        # windowed positions per waveband
        ra = data['ALPHAWIN_J2000_' + WAVEBAND]
        dec = data['DELTAWIN_J2000_' + WAVEBAND]

        delta_ra = (ra - radec_centre[0])*3600.0 * \
                   np.cos(np.deg2rad(radec_centre[1]))
        delta_dec = (dec - radec_centre[1])*3600.0

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

        do_parameter_analysis = False
        if do_parameter_analysis:
            des_parameter_analysis(data=data, index=itest)

        # Measurement image positions per waveband
        ALPHAWIN_J2000 = data['ALPHAWIN_J2000_' + WAVEBAND][itest]
        DELTAWIN_J2000 = data['DELTAWIN_J2000_' + WAVEBAND][itest]

        AWIN_IMAGE = data['AWIN_IMAGE_' + WAVEBAND][itest]
        BWIN_IMAGE = data['BWIN_IMAGE_' + WAVEBAND][itest]
        THETAWIN_IMAGE = data['THETAWIN_IMAGE_' + WAVEBAND][itest]

        X2 = data['X2WIN_IMAGE_' + WAVEBAND][itest]
        Y2 = data['Y2WIN_IMAGE_' + WAVEBAND][itest]
        XY = data['XYWIN_IMAGE_' + WAVEBAND][itest]

        # loop through the sources
        print('Loop through the sources:', len(xdata), WAVEBAND)
        for i, ra, in enumerate(xdata):

            print(i, WAVEBAND, delta_ra[i], delta_dec[i])
            # plot the sources as colored filled circles
            if not plot_ellipses:
                circle = Circle([delta_ra[i], delta_dec[i]],
                                radius,
                                edgecolor='none', facecolor=colors[iband],
                                alpha=alpha)
                plt.gca().add_patch(circle)

            coadd_objects_id = COADD_OBJECTS_ID[i]

            if do_parameter_analysis:
                width = AWIN_IMAGE[i]
                height = BWIN_IMAGE[i]
                angle = THETAWIN_IMAGE[i] + 90.0

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
                    "{:5.1f}".format(XMAX_IMAGE[i] - XMIN_IMAGE[i]),
                    "{:5.1f}".format(YMAX_IMAGE[i] - YMIN_IMAGE[i]),
                    "{:7.1f}".format((XMIN_IMAGE[i] + XMAX_IMAGE[i]) / 2.0),
                    "{:7.1f}".format((YMIN_IMAGE[i] + YMAX_IMAGE[i]) / 2.0),
                    "{:4d}".format(FLAGS[i]))

            # plot as ellipse using the
            width = 1.5 * AWIN_IMAGE[i]
            height = 1.5 * BWIN_IMAGE[i]
            # by eye it looks like angle is reflected
            # maybe need to swap if/when x-axis inverted
            angle = 180.0 - THETAWIN_IMAGE[i]

            print('i, a, b, theta, band, dRA, dDec:',
                  i, height, width, angle, WAVEBAND,
                  delta_ra[i], delta_dec[i])
            print(AWIN_IMAGE[i], BWIN_IMAGE[i], THETAWIN_IMAGE[i])

            a, b, theta = moments2ellipse(X2[i], Y2[i], XY[i])
            print('i, a, b, theta:', i, a, b, theta)


            # https://matplotlib.org/api/_as_gen/matplotlib.patches.Ellipse.html
                  # plot the ellipse marker and edge separately
            ellipse = Ellipse([delta_ra[i], delta_dec[i]],
                width=width/2.0, height=height/2.0, angle=angle,
                edgecolor='none', facecolor=colors[iband],
                alpha=alpha)
            plt.gca().add_patch(ellipse)

            # add solid line to edgecolor for G and Z
            # G and Z chosen since large wavelength range
            # Y not used since lower S/N that Z in DES
            # G is often poorer seeing than Z so can bias size
            if WAVEBAND in ['G', 'Z']:
                edgecolor=colors[iband]
                ellipse = Ellipse([delta_ra[i], delta_dec[i]],
                              fill=False,
                              width=width/2.0, height=height/2.0,
                              angle=angle,
                              facecolor='none',
                              edgecolor=edgecolor,
                              alpha=1.0)

            plt.gca().add_patch(ellipse)


        plt.plot(delta_ra, delta_dec, '.', color=colors[iband],
            label=WAVEBAND + ': ' + str(ndata))

    plt.suptitle('Windowed positions: ' + sourceName)
    plt.xlabel('Delta RA (arc seconds)')
    plt.ylabel('Delta Dec (arc seconds)')
    plt.legend(fontsize='small')
    plotid()

    plotfile = sourceName + '_COADD_radec.png'
    plt.savefig(plotfile)
    #plt.clf()
    print('Saving: ', plotfile)

    if showplot:
        plt.show()

    return plt
