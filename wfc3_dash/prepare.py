"""
Routines to analyze the data from the COSMOS program.

TO DO:

- persistence correction

"""
import astropy
from astropy.io import fits
from astropy.table import Table as table
import drizzlepac
from drizzlepac import astrodrizzle
import os
import glob
import numpy as np

def split_IMA(root='icxe15wwq', PATH='../RAW/'):
    
    """
    Make FLTs from the individual reads.

    root : str
        Root name of the FLT that will be split into individual images.

    PATH : str
        Path to the directory where the clean RAW files are located.
    """
    
    FLAT_F160W = fits.open(os.path.join(os.getenv('iref'),'uc721145i_pfl.fits'))[1].data
    BP_MASK = fits.open('../REF/new_bp_Oct17.fits')[0].data
    
    ima = fits.open(PATH+root+'_ima.fits.gz')
    flt = fits.open(root+'_flt.fits')
    orig_flt = fits.open(PATH+root+'_flt.fits.gz')
    
    NSAMP = ima[0].header['NSAMP']
    sh = ima['SCI',1].shape
    
    cube = np.zeros((NSAMP, sh[0], sh[1]))
    dq = np.zeros((NSAMP, sh[0], sh[1]), dtype=np.int)
        
    time = np.zeros(NSAMP)
    for i in range(NSAMP):
        cube[NSAMP-1-i, :, :] = ima['SCI',i+1].data*ima['TIME',i+1].header['PIXVALUE']
        dq[NSAMP-1-i, :, :] = ima['DQ',i+1].data
        time[NSAMP-1-i] = ima['TIME',i+1].header['PIXVALUE']
    
    diff = np.diff(cube, axis=0)
    dt = np.diff(time)
        
    readnoise_2D = np.zeros((1024,1024))
    readnoise_2D[512: ,0:512] += ima[0].header['READNSEA']
    readnoise_2D[0:512,0:512] += ima[0].header['READNSEB']
    readnoise_2D[0:512, 512:] += ima[0].header['READNSEC']
    readnoise_2D[512: , 512:] += ima[0].header['READNSED']
    readnoise_2D = readnoise_2D**2
        
    for j in range(1, NSAMP-1):
        print('{}_{:02d}_flt.fits'.format(root,j))
        sci = diff[j,:,:]
        exptime = dt[j]
        var = readnoise_2D + sci*FLAT_F160W    
        err = np.sqrt(var)/exptime
        
        flt[0].header['EXPTIME'] = exptime
        flt['SCI'].data = sci[5:-5,5:-5]/exptime
        flt['ERR'].data = err[5:-5,5:-5]
        flt['DQ'].data = dq[j+1][5:-5,5:-5]
        #flt['DQ'].data = dq[-1][5:-5,5:-5]
        flt['DQ'].data[BP_MASK == 1] += 4
        ### trun the 8192 cosmic ray flag to the standard 3096
        flt['DQ'].data[(flt['DQ'].data & 8192) > 0] -= 4096
        ### remove the 32 flag, these are not consistently bad
        flt['DQ'].data[(flt['DQ'].data & 32) > 0] -= 32
        flt['SAMP'].data = np.zeros((1014,1014)) + 1.
        flt['TIME'].data = np.zeros((1014,1014)) + exptime
        flt[0].header['IMA2FLT'] = (1, 'FLT {} extracted from IMA file'.format(j))        

        print('Writing {}_{:02d}_flt.fits'.format(root,j))
        flt.writeto('{}_{:02d}_flt.fits'.format(root,j), clobber=True)
        
def make_pointing_asn(root='icxe15wwq', master_root='icxe15010'):
    
    """
    Makes an new association table for the reads extracted from a given FLT.

    root : str
        Root name of the FLT that was split into individual images. Will be the root name of the new ASN.

    master_root : str
        Master ASN root.
    """
    
    master_asn = fits.open('{}_asn.fits'.format(master_root))
        
    files = glob.glob('{}_*_flt.fits'.format(root))
    nrows = len(files)
    
    #### Primary HDU
    hdu = master_asn[0].copy()
    tbhdu = fits.new_table(master_asn[1].columns, nrows=nrows+1, fill=True)
    for i in range(nrows):
        tbhdu.data[i] = (files[i].split('_flt')[0].upper(), 'EXP-DTH', True)
     
    tbhdu.data[i+1] = (root.upper(), 'PROD-DTH', True)
            
    tbhdu.header = master_asn[1].header.copy()
    tbhdu.header.update('ASN_ID',root)
    tbhdu.header.update('ASN_TAB','{}_asn.fits'.format(root))
    
    #### Create HDUList and write it to output file
    out_fits = fits.HDUList([hdu,tbhdu])
    
    if 'EXTEND' not in hdu.header.keys():
        hdu.header.update('EXTEND', True, after='NAXIS')
                
    print('Writing {}_asn.fits'.format(root))
    out_fits.writeto('{}_asn.fits'.format(root), clobber=True)


def subtract_background_reads(root='icxe15wwq', master_root='icxe15010', subtract=False, reset_stars_dq=False):
    
    """
    Subtract background from the individual reads. Uses the DRZ and SEG images produced in the 
    unicorn FLT background subtraction. 

    root : str
        Root name of the ASN which describes the collection of images.

    master_root : str
        Root of the master ASN.

    subtract : bool
        By default it does not subtract the background but writes it to the header. Set to True to subtract background.

    reser_strs_dq : bool
        If reset_stars_dq = True it will reset cosmic rays within objects to 0 because the centers of stars are flagged.

    """
    
    import stwcs
    import scipy
    import scipy.optimize
        
    asn = utils.ASNFile('{}_asn.fits'.format(root))
    
    #print 'Read files...'
    ref = fits.open('{}_drz_sci.fits'.format(master_root))
    ref_wcs = stwcs.wcsutil.HSTWCS(ref, ext=0)

    seg = fits.open('{}_drz_seg.fits'.format(master_root))    
    seg_data = np.cast[np.float32](seg[0].data)
          
    yi, xi = np.indices((1014,1014))
     
    #### Loop through FLTs
    models = []
    for exp in asn.exposures:
        flt = fits.open('{}_flt.fits'.format(exp)) #, mode='update')
        flt_wcs = stwcs.wcsutil.HSTWCS(flt, ext=1)
        
        if exp == asn.exposures[0]:
            print('Segmentation image: {}_blot.fits'.format(exp))
            blotted_seg = astrodrizzle.ablot.do_blot(seg_data, ref_wcs, flt_wcs, 1, coeffs=True, interp='nearest', sinscl=1.0, stepsize=10, wcsmap=None)         
            
        mask = (blotted_seg == 0) & (flt['DQ'].data == 0) & (flt[1].data > -1) & (xi > 10) & (yi > 10) & (xi < 1004) & (yi < 1004)
        mask &= (flt[1].data < 5*np.median(flt[1].data[mask]))
        data_range = np.percentile(flt[1].data[mask], [2.5, 97.5])
        mask &= (flt[1].data >= data_range[0]) & (flt[1].data <= data_range[1])
        data_range = np.percentile(flt[2].data[mask], [0.05, 99.5])
        mask &= (flt[2].data >= data_range[0]) & (flt[2].data <= data_range[1])
        
        sky_level = np.median(flt[1].data[mask])
        model = flt[1].data*0. + sky_level
        
        # add header keywords of the fit components
        flt = fits.open('{}_flt.fits'.format(exp), mode='update')
        flt[1].header['MDRIZSKY'] =  sky_level 
        if subtract:
            flt[1].data -= model           
            flt[1].header['BG_SUB'] =  'Yes'
        else:
            flt[1].header['BG_SUB'] =  'No'
        
        if reset_stars_dq:
            flagged_stars = ((flt['DQ'].data & 4096) > 0) & (blotted_seg > 0)
            flt['DQ'].data[flagged_stars] -= 4096
        
        flt.flush()
        print('Background subtraction, {}_flt.fits:  {}'.format(exp, sky_level))
        
def fix_cosmic_rays(root='icxe15wwq', master_root = 'icxe15010'):
    
    """ Resets cosmic rays within the seg maps of objects and uses L.A.Cosmic to find them again.

    root : str

    master_root : str
        

    """
    
    from cosmics import cosmics
    import stwcs

    asn = utils.ASNFile('{}_asn.fits'.format(root))
    
    ref = fits.open('{}_drz_sci.fits'.format(master_root))
    ref_wcs = stwcs.wcsutil.HSTWCS(ref, ext=0)

    seg = fits.open('{}_drz_seg.fits'.format(master_root))    
    seg_data = np.cast[np.float32](seg[0].data)
    
    flt_full = fits.open('{}_flt.fits'.format(root))
    flt_full_wcs = stwcs.wcsutil.HSTWCS(flt_full, ext=1)
    
    blotted_seg = astrodrizzle.ablot.do_blot(seg_data, ref_wcs, flt_full_wcs, 1, 
        coeffs=True, interp='nearest', sinscl=1.0, stepsize=10, wcsmap=None)
    
    EXPTIME = flt_full[0].header['EXPTIME']
    SKY = flt_full[0].header['BGCOMP1']    
    
    yi, xi = np.indices((1014,1014))    

    # Build the object :
    #c = cosmics.cosmicsimage(flt_full[1].data*EXPTIME, gain=1.0, readnoise=2.2, 
    #    sigclip = 5.0, sigfrac = 0.4, objlim = 5.0, satlevel=-1., pssl = SKY*EXPTIME)
    c = cosmics.cosmicsimage(flt_full[1].data, gain=1.0, readnoise=0.12, 
        sigclip = 4.0, sigfrac = 0.5, objlim = 10.0, satlevel=-1., pssl = 0.)
    
    # Run the full artillery :
    c.run(maxiter = 4, verbose=True)
    
    for exp in asn.exposures:
        
        flt = fits.open('{}_flt.fits'.format(exp), mode = 'update')
        
        flagged_stars = ((flt['DQ'].data & 4096) > 0) & (blotted_seg > 0)
        flt['DQ'].data[flagged_stars] -= 4096
             
        #new_cr = (c.mask == 1) & ((flt['DQ'].data & 4096) == 0) & ((blotted_seg == 0) | ((blotted_seg > 0) & (flt['SCI'].data < 1.)))  & (xi > 915) & (yi < 295)
        new_cr = (c.mask == 1) & ((flt['DQ'].data & 4096) == 0)
        
        flt['DQ'].data[new_cr] += 4096
        
        flt.flush()
            
        
def align_reads(root='icxe15wwq', threshold=3, final_scale=0.12, refimage='../REF/cosmos-wide_ACS.fits', master_catalog='../REF/IPAC_ACS.fits', align=True, refxcol = 5, refycol = 6):
    
    from drizzlepac import tweakreg

    asn = utils.ASNFile('{}_asn.fits'.format(root))

    catfile = '{}.catfile'.format(root)
    fp = open(catfile,'w')

    #drizzlepac.astrodrizzle.AstroDrizzle('{}_asn.fits'.format(root), output=root, clean=False, context=False, 
    #    preserve=False, skysub=True, driz_separate=True, driz_sep_wcs=True, median=True, blot=True, driz_cr=True, 
    #    driz_cr_corr=True, driz_combine=True)

    for exp in asn.exposures:
        
        flt = fits.open('{}_flt.fits'.format(exp), mode='update')
        if flt[1].header['BG_SUB'] == 'No':
            flt[1].data -= flt[1].header['MDRIZSKY']
            flt[1].header['BG_SUB'] =  'Yes'
            flt.flush()
        
        se = utils.sex.SExtractor()
        se.aXeParams()
        se.copyConvFile()
        se.options['CHECKIMAGE_TYPE'] = 'NONE'
        se.options['FILTER']    = 'N'
        se.options['WEIGHT_IMAGE'] = '{}_flt.fits[1]'.format(exp)
        se.options['WEIGHT_TYPE'] = 'MAP_WEIGHT'
        
        #
        se.params['X_IMAGE'] = True; se.params['Y_IMAGE'] = True
        se.params['MAG_AUTO'] = True
        #
        se.options['CATALOG_NAME'] = '{}_flt.cat'.format(exp)
        se.options['DETECT_THRESH'] = '{}'.format(threshold)
        se.options['ANALYSIS_THRESH'] = '{}' .format(threshold)
        se.options['DETECT_MINAREA'] = '10'
        #
        se.sextractImage('{}_flt.fits[1]'.format(exp))
        utils.sex.sexcatRegions('{}_flt.cat'.format(exp), '{}_flt.reg'.format(exp), format=1)
        
        line = '{0}_flt.fits {0}_flt.cat\n'.format(exp)
        fp.write(line)
        
    fp.close()

    files = glob.glob('{}_*_flt.fits'.format(root))
        
    if align:
        #### Make room for TWEAK wcsname
        for exp in asn.exposures:
            utils.clean_wcsname(flt='{}_flt.fits'.format(exp), wcsname='TWEAK') 
            utils.clean_wcsname(flt='{}_flt.fits'.format(exp), wcsname='OPUS')    
        
        tweakreg.TweakReg(files, refimage=refimage, updatehdr=True, updatewcs=True, catfile=catfile, xcol=2, ycol=3, xyunits='pixels', refcat=master_catalog, refxcol = refxcol, refycol = refycol, refxyunits='degrees', shiftfile=True, fitgeometry='shift',outshifts='{}_shifts.txt'.format(root), outwcs='{}_wcs.fits'.format(root), searchrad=5., tolerance=1., minobj = 5, xoffset = 0.0, yoffset = 0.0, wcsname='TWEAK', interactive=False, residplot='No plot', see2dplot=True, clean=True, headerlet=False, clobber=True)
    
    # AstroDrizzle doesn't like the asn file here: '{}_asn.fits'.format(root)
    # Temporaroly substituting with a list of files
    #files = glob.glob('{}_*_flt.fits'.format(root))
    drizzlepac.astrodrizzle.AstroDrizzle(files, output=root, clean=True, final_scale=final_scale, 
        final_pixfrac=0.8, context=False, resetbits=0, final_bits=576, driz_cr = True, driz_sep_bits=576, 
        preserve=False, wcskey='TWEAK', driz_cr_snr='8.0 5.0', driz_cr_scale = '2.5 0.7', skyuser = 'MDRIZSKY', 
        final_wcs=True)
        
    for exp in asn.exposures:
        flt = fits.open('{}_flt.fits'.format(exp), mode='update')
        flt[1].data += flt[1].header['MDRIZSKY']
        flt[1].header['BG_SUB'] =  'No'
        flt.flush()
    
    
def prep_FLTs(root='icxe15010', refimage='../REF/cosmos-wide_ACS.fits', REF_CAT='../REF/IPAC_ACS.fits'):
    
    
    outshifts = 'shifts_{}.txt'.format(root)
    outwcs = 'shifts_{}_wcs.fits'.format(root)
    
    #drizzlepac.astrodrizzle.AstroDrizzle(root+'_asn.fits', clean=False, context=False, preserve=False, skysub=True, driz_separate=True, driz_sep_wcs=True, median=True, blot=True, driz_cr=True, driz_cr_corr=True, driz_combine=True)
    
    utils.subtract_flt_background(root=root)
    
    drizzlepac.tweakreg.TweakReg(root+'_asn.fits', refimage=refimage, updatehdr = True, updatewcs = True, 
        writecat = False, clean = True, verbose = True, runfile = 'tweakreg.log', 
        wcsname = 'TWEAK', headerlet = False, shiftfile = True, outshifts = outshifts, outwcs = outwcs, 
        refcat = REF_CAT, refxcol = 5, refycol = 6, refxyunits = 'degrees', minobj = 5, searchrad = 1000.0, 
        searchunits = 'pixels', 
        use2dhist = True, see2dplot = False, separation = 0.5, tolerance = 1.0, xoffset = 0.0, yoffset = 0.0, 
        fitgeometry = 'shift', interactive=False, nclip = 3, sigma = 3.0, clobber=True) 
    
    drizzlepac.astrodrizzle.AstroDrizzle(root+'_asn.fits', clean=False, final_pixfrac=1.0, context=False, final_bits=576, resetbits=0, preserve=False, driz_cr_snr='8.0 5.0', driz_cr_scale = '2.5 0.7', wcskey= 'TWEAK')
        
def run_sextractor(mosaic='test1_drz_sci.fits', weight='test1_drz_wht.fits'):
    
    import os
    
    catalog = mosaic.replace('.fits','.cat')
    segmentaion = mosaic.replace('.fits','_seg.fits')
    
    sextr = "sex %s -c gyro.config -CATALOG_NAME %s -MAG_ZEROPOINT %f -BACK_TYPE AUTO,AUTO -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_GAIN Y,Y -WEIGHT_IMAGE %s -GAIN_KEY EXPTIME*CCDGAIN -CHECKIMAGE_NAME %s,%s,%s" %(mosaic, catalog, 25.956, weight, mosaic.replace('.fits','_seg.fits'), mosaic.replace('.fits','_bg.fits'), mosaic.replace('.fits','_sub.fits'))

    os.system(sextr)
    
    
def run_orbit(master_root='icxe15010', RAW_PATH = '../RAW/'):
    
    asn = utils.ASNFile('{}_asn.fits'.format(master_root))
    
    for root in asn.exposures:
        os.system('rsync -av {}/{}_flt.fits.gz .'.format(RAW_PATH, root))
        os.system('gunzip -f {}_flt.fits.gz'.format(root))
        
    prep_FLTs(root=master_root)
    
    for root in asn.exposures:
        split_IMA(root=root)
        make_pointing_asn(root=root, master_root=master_root)
        subtract_background_reads(root=root, master_root=master_root)
        
    fix_cosmic_rays(root=asn.exposures[0], master_root=master_root)
    align_reads(root=asn.exposures[0], align=False)
    for root in asn.exposures[1:]:
        fix_cosmic_rays(root=root, master_root=master_root)
        align_reads(root=root)
    
    #files = glob.glob(master_root[:6]+'*_*_flt.fits')
    #drizzlepac.astrodrizzle.AstroDrizzle(files, output='test3', clean=True, final_scale=0.1, final_pixfrac=0.8, resetbits=0, context=False, preserve=False, skysub = True, skywidth = 0., skystat = '', skylower = None, skyupper = None, skyclip = 0, skylsigma = 0.0, skyusigma = 0.0, skyuser = 'MDRIZSKY', skyfile = '', wcskey = 'TWEAK', driz_separate = False, driz_sep_wcs = False, median = False, blot = False, driz_cr = False, driz_combine = True, final_wht_type = 'IVM', final_kernel = 'square', final_wt_scl = 'exptime', final_fillval = 0,final_bits = 576, final_units = 'cps', final_wcs = True, driz_sep_bits = 576, final_rot=0, final_ra=1.501375000000E+02, final_dec=2.597027777778E+00,driz_cr_snr='8.0 5.0', driz_cr_scale = '2.5 0.7', final_outnx=9100, final_outny=10200)
    
def make_cutouts(root='test6_drz',catalog='test6_drz_sci.cat', DIR='cosmos_wide_stamps/'):
    
    import my_python.mk_region_file
    
    img = fits.open(root+'_sci.fits')
    wht = fits.open(root+'_wht.fits')
    seg = fits.open(root+'_sci_seg.fits')
     
    cat = table.read(catalog, format='ascii.sextractor')
    index = np.where((cat['X_WORLD'] < 150.20664) & (cat['Y_WORLD'] < 2.5588802) & (cat['FLUX_RADIUS'] > 2.))[0]

    my_python.mk_region_file.mk_region_file_from_lists(cat['X_WORLD'][index],cat['Y_WORLD'][index],outfile = 'size', printids='no', color='cyan')
    sz=40.
    
    im_shape = np.shape(img[0].data)

    if not os.path.exists(DIR):
        os.system('mkdir {}'.format(DIR))

    for ii in index:
        print('Making stamp for {}'.format(cat['NUMBER'][ii]))
        x = np.round(cat['X_IMAGE'][ii])
        y = np.round(cat['Y_IMAGE'][ii])
        if (x < sz) or (y < sz) or (x > im_shape[1]-sz) or (y > im_shape[0]-sz):
            continue
        stamp_img = img[0].data[(y-sz):(y+sz),(x-sz):(x+sz)]
        stamp_wht = wht[0].data[(y-sz):(y+sz),(x-sz):(x+sz)]
        stamp_seg = seg[0].data[(y-sz):(y+sz),(x-sz):(x+sz)]

        out_img = fits.PrimaryHDU(stamp_img)
        out_img.writeto('{}cosmos_wide_{:05d}_img.fits'.format(DIR, cat['NUMBER'][ii]))
        out_wht = fits.PrimaryHDU(stamp_wht)
        out_wht.writeto('{}cosmos_wide_{:05d}_wht.fits'.format(DIR, cat['NUMBER'][ii]))
        out_seg = fits.PrimaryHDU(stamp_seg)
        out_seg.writeto('{}cosmos_wide_{:05d}_seg.fits'.format(DIR, cat['NUMBER'][ii]))

    