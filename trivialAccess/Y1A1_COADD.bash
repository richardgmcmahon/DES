
# note the different grammar for id
# 
# coadd_objects_id
# catalogid_g
# imageid_g
#

set -x

stime=$(date '+%s')
echo "Starting program at:  $(date)"

# loop through some tables

MAG_LIMIT="21.0"

#for SN_FIELD in SNE SNS SNX SNC
#for SN_FIELD in COSMOS

RELEASE="Y1A1"

# Y1A1_COADD_OBJECTS_S82
for FIELD in STRIPE82

do

  # create lower case string from SN_FIELD
  sn_field=$( echo $SN_FIELD | tr '[:upper:]' '[:lower:]')
  #sn_field=$( echo $SN_FIELD | tr '[:upper:]' '[:lower:]')

  # note double quotes and lack spaces before or after =
  TABLE="Y1A1_COADD_"${FIELD}

  echo 'TABLE: ' ${TABLE}

  QUERY_ignore="\
    SELECT \
      * \
    FROM \
      ${TABLE}"

   # WHERE \
   #   mag_psf_i < ${MAG_LIMIT}"

  QUERY1="SELECT ra, dec, tilename, run, coadd_objects_id, mag_psf_g, mag_psf_r, mag_psf_i, mag_psf_z, mag_psf_y, magerr_psf_g, magerr_psf_r, magerr_psf_i, magerr_psf_z, magerr_psf_y, mag_model_g, mag_model_r, mag_model_i, mag_model_z, mag_model_y,spread_model_g, spread_model_r, spread_model_i, spread_model_z, spread_model_y, spreaderr_model_g, spreaderr_model_r, spreaderr_model_i, spreaderr_model_z, spreaderr_model_y, flags_g, flags_r, flags_i, flags_z, flags_y, flags_weight_g,flags_weight_r,flags_weight_i,flags_weight_z,flags_weight_y,flags_model_g,flags_model_r,flags_model_i,flags_model_z,flags_model_y, xwin_image_g, ywin_image_g, xwin_image_r, ywin_image_r, xwin_image_i, ywin_image_i, xwin_image_z, ywin_image_z FROM ${TABLE} WHERE mag_psf_z < ${MAG_LIMIT} "


  QUERY2="SELECT ra, dec, tilename, run, coadd_objects_id, mag_psf_g, mag_psf_r, mag_psf_i, mag_psf_z, mag_psf_y, magerr_psf_g, magerr_psf_r, magerr_psf_i, magerr_psf_z, magerr_psf_y, mag_model_g, mag_model_r, mag_model_i, mag_model_z, mag_model_y,spread_model_g, spread_model_r, spread_model_i, spread_model_z, spread_model_y, spreaderr_model_g, spreaderr_model_r, spreaderr_model_i, spreaderr_model_z, spreaderr_model_y, flags_g, flags_r, flags_i, flags_z, flags_y, flags_weight_g,flags_weight_r,flags_weight_i,flags_weight_z,flags_weight_y,flags_model_g,flags_model_r,flags_model_i,flags_model_z,flags_model_y, xwin_image_g, ywin_image_g, xwin_image_r, ywin_image_r, xwin_image_i, ywin_image_i, xwin_image_z, ywin_image_z FROM ${TABLE}"

  QUERY3="\
    SELECT \
      * \
    FROM \
      ${TABLE}"

  #QUERY4="SELECT ra, dec, tilename, run, coadd_objects_id, mag_psf_g, mag_psf_r, mag_psf_i, mag_psf_z, mag_psf_y FROM ${TABLE}"

  QUERY5="SELECT ra, dec, tilename, run, coadd_objects_id, mag_psf_g, mag_psf_r, mag_psf_i, mag_psf_z, mag_psf_y FROM ${TABLE} WHERE mag_psf_z < ${MAG_LIMIT} "
  QUERY=$QUERY1

  echo "QUERY: " ${QUERY}

  #read -p "Press [Any] key to continue..."

  #FITSFILE="/data/des/SVA1/${SN_FIELD}/coadd_${sn_field}_ilt21p0_tmp.fits"
  #FITSFILE="/data/des/SVA1/${SN_FIELD}/coadd_${sn_field}_ilt21p0.fits"

  FITSFILE="/data/desardata/Y1A1/${FIELD}/Y1A1_COADD_${FIELD}_thin_ilt21.fits"
  FITSFILE="/data/desardata/Y1A1/${FIELD}/Y1A1_COADD_${FIELD}_thin_zlt21.fits"

  #FITSFILE="/data/desardata/SVA1/${SN_FIELD}/sva1_coadd_${sn_field}_thin.fits"

  echo $FITSFILE

  # note the spaces after the protected double quotes
  trivialAccess --command="write_fits --file=$FITSFILE \" $QUERY \" "

  now=$(date '+%s')
  etime=$((now - stime))
  echo "Elapsed time (secs) : " $etime

done


now=$(date '+%s')
etime=$((now - stime))
echo "Elapsed time (secs) : " $etime
