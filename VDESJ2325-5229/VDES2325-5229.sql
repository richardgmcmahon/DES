SELECT 
    o.*,
    i.fwhm, i.psf_fwhm, i.psf_beta, 
    i.zeropoint, i.sigma_zeropoint,
    e.mjd_obs, e.nite, 
    c.id AS cid, 
    c.ccd, c.tilename AS ctilename, c.objects, 
    c.path AS catalogpath, 
    i.id AS iid, 
    i.path AS imagepath 
FROM
    Y1A1_FINALCUT  o 
    /* Y1A1_OBJECTS  o */
    /* Y1A1_COADD_OBJECTS  o */
    INNER JOIN Y1A1_CATALOG  c 
        ON  o.catalogid = c.id
    INNER JOIN Y1A1_IMAGE   i
        ON  o.imageid = i.id
    INNER JOIN EXPOSURE  e
        ON  c.exposureid = e.id
WHERE
    (o.ra BETWEEN 351.40784 AND 351.43567) 
    AND (o.dec BETWEEN -52.49601 AND -52.47906);
