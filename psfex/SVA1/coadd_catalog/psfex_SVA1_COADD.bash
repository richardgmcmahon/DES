# create catalog Sextractor file with vignettes and then run PSFex

# this will be a module in future
. /opt/ioa/setup/setup_heasoft.bash 
#fcopy -v
#funpack -H
#imcopy -v



. /opt/ioa/Modules/default/init/bash
module purge
module load sextractor/2.19.5
module load psfex/3.17.1

module list

RUN_SEXTRACTOR=true
RUN_PSFEX=true

RELEASE=SVA1
RELEASE=Y1A1

for WAVEBAND in g i
#for WAVEBAND in g r i z Y
#for WAVEBAND in r z Y

do

#for TILE in DES2327-5248 
#for TILE in DES0449-4706 
for TILE in DES0406-5414 
#DES0449-4748 
#for TILE in DES0408-5353 DES2327-5248 

do

DATAPATH_IN=/data/desardata/PSFEX/${RELEASE}/${TILE}/
DATAPATH_OUT=/data/desardata/PSFEX/${RELEASE}/${TILE}/R1

mkdir ${DATAPATH_OUT}

if [ -d ${DATAPATH_OUT} ]
then
    echo  ${DATAPATH_OUT} "does not exist; creating " ${DATAPATH_OUT}
    mkdir ${DATAPATH_OUT}
    echo 
fi

INFILE=${DATAPATH_IN}/${TILE}_${WAVEBAND}.fits
IMAGE=${INFILE}[1]
WEIGHTMAP=${INFILE}[2]

# check is a .fits file or .fits.fz file exists
if [ -f $INFILE ];
then
   echo "File $FILE exists."
else
   echo "File $FILE does not exist."
   if [ -f $INFILE.fz ];
   echo "File $FILE exists."
   then
      echo "File $FILE.fz exists."
      #funpack $INFILE.fz 
      imcopy $INFILE.fz  $INFILE
   else
      echo "File $FILE.fz does not exist."
   fi
fi

OUTFILE=${DATAPATH_OUT}/${TILE}_${WAVEBAND}_psfcat.fits

echo 'Input image: ' ${INFILE} 
echo

time0=$(date '+%s')


if $RUN_SEXTRACTOR; then

echo "Starting Sextractor at:  $(date)"


# just to see the version info
sex --version

now=$(date '+%s')
echo "Elapsed time (secs) : " $((now - time0))
echo 'Sextractor starting'

sex ${IMAGE} -c default.sex -CATALOG_NAME ${OUTFILE} -CATALOG_TYPE FITS_LDAC -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE ${WEIGHTMAP} -PARAMETERS_NAME sex.param_psfex  -FILTER_NAME sex.conv -STARNNW_NAME sex.nnw  -SATUR_LEVEL 65000 -DETECT_MINAREA 3 

echo 'Sextractor completed'
now=$(date '+%s')
etime=$((now - time0))
echo "Elapsed time (secs) : " $etime

fi

if $RUN_PSFEX; then

echo "Starting psfex at:  $(date)"

INFILE_PSFEX=${DATAPATH_OUT}/${TILE}_${WAVEBAND}_psfcat.fits
OUTCAT_PSFEX=${DATAPATH_OUT}/${TILE}_${WAVEBAND}_psfcat-starlist.fits

psfex ${INFILE_PSFEX} \
  -c default_20160202.psfex \
  -WRITE_XML Y \
  -OUTCAT_NAME ${OUTCAT_PSFEX} 

#  -CHECKPLOT_NAME ${DATAPATH_OUT}/fwhm, ${DATAPATH_OUT}/ellipticity, ${DATAPATH_OUT}/counts, ${DATAPATH_OUT}/countfrac, ${DATAPATH_OUT}/chi2, ${DATAPATH_OUT}/residuals

mv *${TILE}*png ${DATAPATH_OUT}/.

echo 'PSFEx completed'
now=$(date '+%s')
etime=$((now - time0))
echo "Elapsed time (secs) : " $etime

fi

done

done