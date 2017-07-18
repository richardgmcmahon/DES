from __future__ import print_function, division


import numpy as np

from astropy.table import Table, join, hstack, vstack
from astropy.io import ascii

def mad(data, median=None, sigma=False):
    if median is None:
        median=np.median(data)
        mad=np.median (abs(data-median))
    if sigma:
        mad=mad/0.6745
    return mad



def main(datapath=None, filename=None, release=None, tile=None):

    t=Table.read('/data/desardata/'+release+'/'+tile+'/'+tile+'_g_cat.fits')
    columns=t.columns
    t=Table.read('/data/desardata/'+release+'/'+tile+'/'+tile+'_merged_cat.fits')
    n=0
    number=[]
    name=[]
    samesies=[]
    gmin=[]
    gmax=[]
    gmedian=[]
    gmad=[]
    rmin=[]
    rmax=[]
    rmedian=[]
    rmad=[]
    imin=[]
    imax=[]
    imedian=[]
    imad=[]
    zmin=[]
    zmax=[]
    zmedian=[]
    zmad=[]
    Ymin=[]
    Ymax=[]
    Ymedian=[]
    Ymad=[]
    for column in columns:
        if column=='FLUX_APER' or column=='FLUXERR_APER' or column=='MAG_APER' or column=='MAGERR_APER':
            col=column
            for band in bands:
                print('band:', band)
                a1=[]
                a2=[]
                a3=[]
                a4=[]
                a5=[]
                a6=[]
                a7=[]
                a8=[]
                a9=[]
                a10=[]
                a11=[]
                a12=[]
                for o in xrange(len(t)):
                    a1.append(t[column+band][o][0])
                    a2.append(t[column+band][o][1])
                    a3.append(t[column+band][o][2])
                    a4.append(t[column+band][o][3])
                    a5.append(t[column+band][o][4])
                    a6.append(t[column+band][o][5])
                    a7.append(t[column+band][o][6])
                    a8.append(t[column+band][o][7])
                    a9.append(t[column+band][o][8])
                    a10.append(t[column+band][o][9])
                    a11.append(t[column+band][o][10])
                    a12.append(t[column+band][o][11])

                t[column+'_1'+band]=np.array(a1)
                t[column+'_2'+band]=np.array(a2)
                t[column+'_3'+band]=np.array(a3)
                t[column+'_4'+band]=np.array(a4)
                t[column+'_5'+band]=np.array(a5)
                t[column+'_6'+band]=np.array(a6)
                t[column+'_7'+band]=np.array(a7)
                t[column+'_8'+band]=np.array(a8)
                t[column+'_9'+band]=np.array(a9)
                t[column+'_10'+band]=np.array(a10)
                t[column+'_11'+band]=np.array(a11)
                t[column+'_12'+band]=np.array(a12)

            for c in xrange(12):
                column=col+'_'+str(c+1)
                print('column:', column)
                number.append(n+1)
                name.append(column)
                test1=t[column+'_G']-t[column+'_R']
                test2=t[column+'_G']-t[column+'_I']
                test3=t[column+'_G']-t[column+'_Z']
                test4=t[column+'_G']-t[column+'_Y']

                if min(test1)==0. and min(test2)==0. and min(test3)==0. and min(test4)==0. and max(test1)==0. and max(test2)==0. and max(test3)==0. and max(test4)==0. and np.median(test1)==0. and np.median(test2)==0. and np.median(test3)==0. and np.median(test4)==0.:
                    samesies.append('det')
                else:
                    samesies.append('mes')
                    gmin.append(min(t[column+'_G']))
                    gmax.append(max(t[column+'_G']))
                    gmedian.append(np.median(t[column+'_G']))
                    gmad.append(mad(t[column+'_G']))
                    rmin.append(min(t[column+'_R']))
                    rmax.append(max(t[column+'_R']))
                    rmedian.append(np.median(t[column+'_R']))
                    rmad.append(mad(t[column+'_R']))
                    imin.append(min(t[column+'_I']))
                    imax.append(max(t[column+'_I']))
                    imedian.append(np.median(t[column+'_I']))
                    imad.append(mad(t[column+'_I']))
                    zmin.append(min(t[column+'_Z']))
                    zmax.append(max(t[column+'_Z']))
                    zmedian.append(np.median(t[column+'_Z']))
                    zmad.append(mad(t[column+'_Z']))
                    Ymin.append(min(t[column+'_Y']))
                    Ymax.append(max(t[column+'_Y']))
                    Ymedian.append(np.median(t[column+'_Y']))
                    Ymad.append(mad(t[column+'_Y']))
        else:
            number.append(n+1)
            name.append(column)
            test1=t[column+'_G']-t[column+'_R']
            test2=t[column+'_G']-t[column+'_I']
            test3=t[column+'_G']-t[column+'_Z']
            test4=t[column+'_G']-t[column+'_Y']
            if min(test1)==0. and min(test2)==0. and min(test3)==0. and min(test4)==0. and max(test1)==0. and max(test2)==0. and max(test3)==0. and max(test4)==0. and np.median(test1)==0. and np.median(test2)==0. and np.median(test3)==0. and np.median(test4)==0.:
                samesies.append('det')
            else:
                samesies.append('mes')
            gmin.append(min(t[column+'_G']))
            gmax.append(max(t[column+'_G']))
            gmedian.append(np.median(t[column+'_G']))
            gmad.append(mad(t[column+'_G']))
            rmin.append(min(t[column+'_R']))
            rmax.append(max(t[column+'_R']))
            rmedian.append(np.median(t[column+'_R']))
            rmad.append(mad(t[column+'_R']))
            imin.append(min(t[column+'_I']))
            imax.append(max(t[column+'_I']))
            imedian.append(np.median(t[column+'_I']))
            imad.append(mad(t[column+'_I']))
            zmin.append(min(t[column+'_Z']))
            zmax.append(max(t[column+'_Z']))
            zmedian.append(np.median(t[column+'_Z']))
            zmad.append(mad(t[column+'_Z']))
            Ymin.append(min(t[column+'_Y']))
            Ymax.append(max(t[column+'_Y']))
            Ymedian.append(np.median(t[column+'_Y']))
            Ymad.append(mad(t[column+'_Y']))

        n+=1

    tout=Table([number,name,samesies,gmin,rmin,imin,zmin,Ymin,gmax,rmax,imax,zmax,Ymax,gmedian,rmedian,imedian,zmedian,Ymedian,gmad,rmad,imad,zmad,Ymad],names=['number','column','imagedetect','gmin','rmin','imin','zmin','Ymin','gmax','rmax','imax','zmax','Ymax','gmedian','rmedian','imedian','zmedian','Ymedian','gmad','rmad','imad','zmad','Ymad'])
    tout.write(tile+'_'+release+'_columns.fits',overwrite=True)
    ascii.write(tout,tile+'_'+release+'_columns.txt')
    print('###########################################')


if __name__ == '__main__':


    tile='DES0453-4457'
    releases=['SVA1','Y1A1']
    #releases=['Y1A1','SVA1']
    bands=['_G','_R','_I','_Z','_Y']


    for release in releases:

        main(datapath=None, filename=None, release=release, tile=tile)
