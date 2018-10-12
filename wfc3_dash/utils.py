import os
from astropy.io import fits
from urllib.request import urlretrieve 


def get_flat(file_name):
    '''
    Will check if user has proper reference file directories 
    and files. Will also return flat field file appropriate for 
    the input file. 

    Parameters
    ----------
    file_name : string
        File name of input IMA. 

    Returns
    ----------
    reffile_name : string
        File name of flat field for that file. 

    '''
    os.environ['iref'] = '~/iref/'
    if not os.path.exists('iref'):
        os.mkdir('iref')
    
    base_url = 'https://hst-crds.stsci.edu/unchecked_get/references/hst/'
    
    with fits.open(file_name) as fitsfile:
        reffile_name = fitsfile[0].header['PFLTFILE'].replace('$', '/')
        if not os.path.exists(reffile_name):
            urlretrieve(base_url + os.path.basename(reffile_name), reffile_name)

    return reffile_name

def footprints_plot(root='icxe15010'):
    
    #    import unicorn.survey_paper as sup  #Mario: commenting out the sup dependencies. Will fix this later
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    
    if root == 'icxe15010':
        aspect = 1.75
        xlim = [150.265, 150.157]
        ylim = [2.45, 2.64]
        xticklab = [r'$10^\mathrm{h}01^\mathrm{m}00^\mathrm{s}$', r'$10^\mathrm{h}00^\mathrm{m}45^\mathrm{s}$']
        #xtickv = [sup.degrees(10,01,00, hours=True),sup.degrees(10,00,45, hours=True)]
        yticklab = [r'$+02^\circ30^\prime00^{\prime\prime}$',r'$+02^\circ35^\prime00^{\prime\prime}$']
        #ytickv = [sup.degrees(2, 30, 00, hours=False),sup.degrees(2, 35, 00, hours=False)]
        label = 'COSMOS-15'
        factor=10.

    if root == 'icxe16010':
        aspect=0.9
        xlim = [150.265, 150.1]
        ylim = [2.607, 2.74]
        xticklab = [r'$10^\mathrm{h}01^\mathrm{m}00^\mathrm{s}$', r'$10^\mathrm{h}00^\mathrm{m}45^\mathrm{s}$',r'$10^\mathrm{h}00^\mathrm{m}30^\mathrm{s}$']
        #xtickv = [sup.degrees(10,01,00, hours=True),sup.degrees(10,00,45, hours=True),sup.degrees(10,00,30, hours=True)]
        yticklab = [r'$+02^\circ38^\prime00^{\prime\prime}$',r'$+02^\circ40^\prime00^{\prime\prime}$', r'$+02^\circ42^\prime00^{\prime\prime}$', r'$+02^\circ44^\prime00^{\prime\prime}$']
        #ytickv = [sup.degrees(2, 38, 00, hours=False),sup.degrees(2, 40, 00, hours=False),sup.degrees(2, 42, 00, hours=False),sup.degrees(2, 44, 00, hours=False)]
        label='COSMOS-16'
        factor=20.
    
    if root == 'icxe17010':
        aspect=1.4
        xlim = [150.2, 150.06]
        ylim = [2.52, 2.72]
        xticklab = [r'$10^\mathrm{h}00^\mathrm{m}45^\mathrm{s}$', r'$10^\mathrm{h}00^\mathrm{m}30^\mathrm{s}$',r'$10^\mathrm{h}00^\mathrm{m}15^\mathrm{s}$']
        #xtickv = [sup.degrees(10,00,45, hours=True),sup.degrees(10,00,30, hours=True),sup.degrees(10,00,15, hours=True)]
        yticklab = [r'$+02^\circ35^\prime00^{\prime\prime}$',r'$+02^\circ40^\prime00^{\prime\prime}$']
        #ytickv = [sup.degrees(2, 35, 00, hours=False),sup.degrees(2, 40, 00, hours=False)]
        label='COSMOS-17'
        factor=240.

    if root == 'icxe18010':
        aspect=1.577
        xlim = [150.14, 150.01]
        ylim = [2.53, 2.735]
        xticklab = [r'$10^\mathrm{h}00^\mathrm{m}30^\mathrm{s}$', r'$10^\mathrm{h}00^\mathrm{m}20^\mathrm{s}$',r'$10^\mathrm{h}00^\mathrm{m}10^\mathrm{s}$']
        #xtickv = [sup.degrees(10,00,30, hours=True),sup.degrees(10,00,20, hours=True),sup.degrees(10,00,10, hours=True)]
        yticklab = [r'$+02^\circ35^\prime00^{\prime\prime}$',r'$+02^\circ40^\prime00^{\prime\prime}$']
        #ytickv = [sup.degrees(2, 35, 00, hours=False),sup.degrees(2, 40, 00, hours=False)]
        label='COSMOS-18'
        factor=240.
    
    
    
        #   fig = unicorn.catalogs.plot_init(square=True, xs=5., aspect=aspect,  # MArio: changed this to a regular matplotlib call (below), until we fix the unicorn module dependency
        #fontsize=8, left=0.18, right=0.02, bottom=0.10, top=0.10)

    fig = plt.figure(square=True, xs=5., aspect=aspect, 
        fontsize=8, left=0.18, right=0.02, bottom=0.10, top=0.10)
    ax = fig.add_subplot(111)
    jet = cm = plt.get_cmap('jet')
    cNorm = colors.Normalize(vmin=0, vmax=9)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)    
    
    reg_file = root+'_asn.reg'
    
    poly = []
    with open(reg_file) as f:
        for line in f:
            if not line.startswith('fk5'):
                region = line.split('#')[0]
                #            poly.append(sup.polysplit(region=region, get_shapely=True))

    shifts = table.read('shifts_{}.txt'.format(root), format='ascii', 
        names=('file','x','y','rot','scale','x_rms','y_rms'))
        
    cc = 0
    xcen_all = []
    ycen_all = []
    for j,(pp, x_off, y_off, file) in enumerate(zip(poly, shifts['x'], shifts['y'], shifts['file'])):
        cc += 1.
        color = scalarMap.to_rgba(cc)
        x, y = pp.exterior.xy
        flt = fits.open(file)
        xcen = flt[1].header['CRVAL1O']
        ycen = flt[1].header['CRVAL2O']
        x_off = (flt[1].header['CRVAL1B']-flt[1].header['CRVAL1O'])*20.
        y_off = (flt[1].header['CRVAL2B']-flt[1].header['CRVAL2O'])*20.
        #print file, xcen, xcen+x_off, ycen, ycen+y_off
        #xcen = (np.mean(x[:-1]))
        #ycen = (np.mean(y[:-1]))
        xcen_all.append(xcen)
        ycen_all.append(ycen)
        ax.plot(x,y,'-', color=color)
        #ax.annotate("",xy=(xcen+(x_off*0.12)/factor, ycen+(y_off*0.12)/factor), xytext=(xcen, ycen), 
        #    arrowprops=dict(arrowstyle='->', color=color))
        #ax.plot([xcen, xcen+x_off], [ycen, ycen+y_off], '-')
        ax.annotate("",xy=(xcen+x_off, ycen+y_off), xytext=(xcen, ycen), 
            arrowprops=dict(arrowstyle='->', color=color))

    ax.plot(xcen_all, ycen_all, '+:', markersize=10., color='0.5', alpha=0.5) 
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)     
    ax.set_xticklabels(xticklab)
    xtick = ax.set_xticks(xtickv)
    ax.set_yticklabels(yticklab)
    ytick = ax.set_yticks(ytickv)
    ax.set_title(label)       
    plt.show(block=False)
    
    fig.savefig('footprint_{}.png'.format(label.lower()), dpi=200, transparent=False)          


class ASNFile(object):
    """
    ASNFile()
    
        Class for handling ASN fits files.
    
        >>> asn = ASNFile(file='ib3701050_asn.fits', grow=False)
        >>> asn.exposures
        ['ib3701ryq', 'ib3701sbq', 'ib3701sdq', 'ib3701sqq']
        >>> asn.product
        'IB3701050'

        If grow=True, allow file rootnames to be 20 characters rather than 14.
     """
    def _read_asn_file(self, grow=True):
        """
        _read_asn_file(self)
        
            Read an ASN FITS file (self.file).
        """
        import numpy as np
        from warnings import warn
        
        self.in_fits = pyfits.open(self.file)
        
        data = self.in_fits[1].data
        
        if grow:
            #### Allow more characters in the MEMNAME column
            memname = pyfits.Column(name='MEMNAME', format='40A', array=self.in_fits[1].columns[0].array.astype('S40'), disp='A40')
            memtype = self.in_fits[1].columns[1]
            memprsnt = self.in_fits[1].columns[2]
            coldefs = pyfits.ColDefs([memname, memtype, memprsnt])
            try:
                #print 'from_columns'
                hdu = pyfits.BinTableHDU.from_columns(coldefs)
            except:
                #print 'fail pyfits'
                hdu = pyfits.new_table(coldefs)
                
            hdu.header = self.in_fits[1].header
            hdu.header['TFORM1'] = '40A'
            hdu.header['TDISP1'] = 'A40'
            hdu.header['NAXIS1'] += 26
            self.in_fits[1] = hdu        
    
        data = self.in_fits[1].data
        #print data
        
        self.header = self.in_fits[0].header
        
        names = data.field('MEMNAME')
        types = data.field('MEMTYPE')
        
        ##### Exposures
        #exp_idx  = np.where(types == 'EXP-DTH')
        exp_idx = types == 'EXP-DTH'
        #### check if MEMTYPE starts with EXP, have other cases where type is "EXP-RP#"
        for ii, type in enumerate(types):
            if types[ii].startswith('EXP'):
                exp_idx[ii] = True
                
        if exp_idx.sum() == 0:
            warn('ASN file %s has no EXP-DTH items')
        else:
            self.exposures = []
            for exp in names[exp_idx]:
                self.exposures.append(exp.lower())
        
        ##### Products
        prod_idx = np.where(types == 'PROD-DTH')
        if prod_idx[0].shape[0] != 1:
            warn('ASN file %s has N != 1 PROD-DTH items' %self.file )
            self.product = None
        else:
            self.product = names[prod_idx[0]][0].upper()
    
    
    def __init__(self, file=None, grow=True):
        self.file = file
        self.exposures = []
        self.product = None
        if file:
            self._read_asn_file(grow=grow)
    
    
    def write(self, out_file=None, clobber=True):
        """
        write(self,out_file=None, clobber=True)

            out_file='self' writes to `self.file`.
    
        """
        if not out_file:
            print("USAGE:: asn.write(out_file='output_asn.fits')")
        else:
            if out_file == 'self':
                out_file = self.file
            
            nexp  = self.exposures.__len__()
            if self.product:
                nprod = 1
            else:
                nprod = 0
            nrows = nexp + nprod
            #### Primary HDU
            hdu = self.in_fits[0].copy()
            #### BinTable HDU
            tbhdu = pyfits.new_table(self.in_fits[1].columns, nrows=nrows, fill=True)
            for i in range(nexp):
                tbhdu.data[i] = (self.exposures[i].upper(), 'EXP-DTH', True)
            if nprod > 0:
                tbhdu.data[i+1] = (self.product, 'PROD-DTH', True)
            
            tbhdu.header = self.in_fits[1].header.copy()
            tbhdu.header.update('ASN_ID',out_file.split('_asn.fits')[0])
            tbhdu.header.update('ASN_TAB',out_file)
            #### Create HDUList and write it to output file
            self.out_fits = pyfits.HDUList([hdu,tbhdu])
            if 'EXTEND' not in list(hdu.header.keys()):
                hdu.header.update('EXTEND', True, after='NAXIS')
                
            self.out_fits.writeto(out_file, clobber=clobber)
    
    def showContents(self):
        """
        showContents()
        
            >>> x = ASNFile(file='ib3702060_asn.fits')
            >>> x.showContents()
            1   ib3703uxq    EXP-DTH      yes
            2   ib3703vaq    EXP-DTH      yes
            3   ib3703vcq    EXP-DTH      yes
            4   ib3703vsq    EXP-DTH      yes
            5   IB3703050    PROD-DTH     yes
        """
        if self.exposures.__len__() > 0:
            for i,exp in enumerate(self.exposures):
                print('%5d   %s    EXP-DTH      yes' %(i+1,exp))
            print('%5d   %s    PROD-DTH     yes' %(i+2,self.product))
    
    def append(self, new):
        """
        append(self, new)
        
            `new` must be an instance of ASNFile.
        
            `new.exposures` are added to the `self.exposures` list.
        """
        from warnings import warn
        if not isinstance(new,self.__class__):
            warn("argument is not an instance of ASNFile()")
        else:
            self.exposures.extend(new.exposures)
    
def clean_wcsname(flt='ibhj15wyq_flt.fits', wcsname='TWEAK', ACS=False, WFPC2=False):
    """
    Workaround for annoying TweakReg feature of not overwriting WCS solns
    """
    im = pyfits.open(flt, mode='update')
    if ACS:
        exts = [1,4]
    elif WFPC2:
        exts = [1,2,3,4]
    else:
        exts = [1]
    
    for ext in exts:
        header = im[ext].header
        for key in header:
            if key.startswith('WCSNAME'):
                if header[key] == wcsname:
                    wcs_ext = key[-1]
                    if key == 'WCSNAME':
                        header[key] = 'X' + wcsname+'X'
                        #im.flush()
        #
        for key in ['WCSNAME', 'WCSAXES', 'CRPIX1', 'CRPIX2', 'CDELT1', 'CDELT2', 'CUNIT1', 'CUNIT2', 'CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2', 'LONPOLE', 'LATPOLE', 'CRDER1', 'CRDER2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', 'FITNAME', 'NMATCH', 'RMS_RA', 'RMS_DEC']:
            try:
                header.remove(key+wcs_ext)
            except:
                #print key
                pass
    
    im.flush()

def subtract_flt_background(root='GOODN-N1-VBA-F105W', scattered_light=False, sex_background=False, order=2):
    """
    Subtract polynomial background
    """
    import scipy.optimize
    
    import astropy.units as u
    
    from astropy.table import Table as table
    
    import stwcs
    from stwcs import updatewcs
    
    import drizzlepac
    from drizzlepac import astrodrizzle, tweakreg, tweakback
    
    import threedhst
    
    asn = threedhst.utils.ASNFile(root+'_asn.fits')
    for exp in asn.exposures:
        updatewcs.updatewcs('%s_%s.fits' %(exp, 'flt'))

    if not os.path.exists('%s_drz_sci.fits' %(root)):        
        if len(asn.exposures) == 1:
            drizzlepac.astrodrizzle.AstroDrizzle(root+'_asn.fits', clean=False, context=False, preserve=False, skysub=True, driz_separate=False, driz_sep_wcs=False, median=False, blot=False, driz_cr=False, driz_cr_corr=False, driz_combine=True)
        else:
            drizzlepac.astrodrizzle.AstroDrizzle(root+'_asn.fits', clean=False, context=False, preserve=False, skysub=True, driz_separate=True, driz_sep_wcs=True, median=True, blot=True, driz_cr=True, driz_cr_corr=True, driz_combine=True)
    
    se = threedhst.sex.SExtractor()
    se.options['WEIGHT_IMAGE'] = '%s_drz_wht.fits' %(root)
    se.options['WEIGHT_TYPE'] = 'MAP_WEIGHT'
    se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION,BACKGROUND'
    se.options['CHECKIMAGE_NAME'] = '%s_drz_seg.fits,%s_drz_bkg.fits' %(root, root)
    se.options['BACK_TYPE'] = 'AUTO'
    se.options['BACK_SIZE'] = '256'
    #
    se.params['X_IMAGE'] = True; se.params['Y_IMAGE'] = True
    se.params['MAG_AUTO'] = True
    #
    se.options['CATALOG_NAME'] = '%s_drz_sci.cat' %(root)
    se.options['FILTER'] = 'Y'
    se.copyConvFile()
    se.options['FILTER_NAME'] = 'gauss_4.0_7x7.conv'
    se.options['DETECT_THRESH'] = '0.8'
    se.options['ANALYSIS_THRESH'] = '0.8'
    #
    se.options['MEMORY_OBJSTACK'] = '30000'
    se.options['MEMORY_PIXSTACK'] = '3000000'
    se.options['MEMORY_BUFSIZE'] = '2048'
    
    se.sextractImage('%s_drz_sci.fits' %(root))
    #threedhst.sex.sexcatRegions('%s_flt.cat' %(exp), '%s_flt.reg' %(exp), format=1)
    
    #### Blot segmentation map to FLT images for object mask
    asn = threedhst.utils.ASNFile('%s_asn.fits' %(root))
    
    #print 'Read files...'
    ref = pyfits.open('%s_drz_sci.fits' %(root))
    ref_wcs = stwcs.wcsutil.HSTWCS(ref, ext=0)

    seg = pyfits.open('%s_drz_seg.fits' %(root))    
    #### Fill ref[0].data with zeros for seg mask
    #seg_data = ref[0].data
    #seg_data[seg[0].data == 0] = 0
    seg_data = np.cast[np.float32](seg[0].data)
    
    bkg_data = pyfits.open('%s_drz_bkg.fits' %(root))[0].data
      
    yi, xi = np.indices((1014,1014))
    if scattered_light:        
        bg_components = np.ones((4,1014,1014))
        bg_components[1,:,:] = xi/1014.*2
        bg_components[2,:,:] = yi/1014.*2
        bg_components[3,:,:] = pyfits.open(os.getenv('THREEDHST') + '/CONF/G141_scattered_light.fits')[0].data
        #### Use flat-field itself for images affected by full-field 
        #### persistence from the tungsten lamp
        if scattered_light == 2:
            bg_components[3,:,:] = pyfits.open(os.getenv('iref') + 'flat_UDF_F140W_v0.fits')[1].data[5:-5,5:-5]
            
        NCOMP=4
    else:
        # bg_components = np.ones((3,1014,1014))
        # bg_components[1,:,:] = xi/1014.*2
        # bg_components[2,:,:] = yi/1014.*2
        # NCOMP=3
        #
        if order == 2:
            NCOMP=6
            bg_components = np.ones((NCOMP,1014,1014))
            bg_components[1,:,:] = (xi-507)/507.
            bg_components[2,:,:] = (yi-507)/507.
            bg_components[3,:,:] = ((xi-507)/507.)**2
            bg_components[4,:,:] = ((yi-507)/507.)**2
            bg_components[5,:,:] = (xi-507)*(yi-507)/507.**2
        else:
            NCOMP=3
            bg_components = np.ones((NCOMP,1014,1014))
            bg_components[1,:,:] = (xi-507)/507.
            bg_components[2,:,:] = (yi-507)/507.
            
    bg_flat = bg_components.reshape((NCOMP,1014**2))
    
    #### Loop through FLTs, blotting reference and segmentation
    models = []
    for exp in asn.exposures:
        flt = pyfits.open('%s_flt.fits' %(exp)) #, mode='update')
        flt_wcs = stwcs.wcsutil.HSTWCS(flt, ext=1)
        
        ### segmentation        
        print('Segmentation image: %s_blot.fits' %(exp))
        blotted_seg = astrodrizzle.ablot.do_blot(seg_data+0, ref_wcs, flt_wcs, 1, coeffs=True, interp='nearest', sinscl=1.0, stepsize=10, wcsmap=None)
        
        blotted_bkg = 0.
        if sex_background:
            blotted_bkg = astrodrizzle.ablot.do_blot(bkg_data+0, ref_wcs, flt_wcs, 1, coeffs=True, interp='nearest', sinscl=1.0, stepsize=10, wcsmap=None)
            flt[1].data -= blotted_bkg
            
        mask = (blotted_seg == 0) & (flt['DQ'].data == 0) & (flt[1].data > -1) & (xi > 10) & (yi > 10) & (xi < 1004) & (yi < 1004) 
        mask &= np.isfinite(flt[1].data) & np.isfinite(flt[2].data)
        mask &= (flt[1].data < 5*np.median(flt[1].data[mask]))
        data_range = np.percentile(flt[1].data[mask], [2.5, 97.5])
        mask &= (flt[1].data >= data_range[0]) & (flt[1].data <= data_range[1])
        data_range = np.percentile(flt[2].data[mask], [0.5, 99.5])
        mask &= (flt[2].data >= data_range[0]) & (flt[2].data <= data_range[1])
        
        ### Least-sq fit for component normalizations
        data = flt[1].data[mask].flatten()
        wht = (1./flt[2].data[mask].flatten())**2
        templates = bg_flat[:, mask.flatten()]
        p0 = np.zeros(NCOMP)
        p0[0] = np.median(data)
        obj_fun = threedhst.grism_sky.obj_lstsq
        print('XXX: %d' %(mask.sum()))
        popt = scipy.optimize.leastsq(obj_fun, p0, args=(data, templates, wht), full_output=True, ftol=1.49e-8/1000., xtol=1.49e-8/1000.)
        xcoeff = popt[0]
        model = np.dot(xcoeff, bg_flat).reshape((1014,1014))
        models.append(model)
        
        # add header keywords of the fit components
        flt = pyfits.open('%s_flt.fits' %(exp), mode='update')
        flt[1].data -= model+blotted_bkg
        for i in range(NCOMP):
            if 'BGCOMP%d' %(i+1) in flt[0].header:
                flt[0].header['BGCOMP%d' %(i+1)] += xcoeff[i]
            else:
                flt[0].header['BGCOMP%d' %(i+1)] = xcoeff[i]                
        
        flt.flush()
        coeff_str = '  '.join(['%.4f' %c for c in xcoeff])
        threedhst.showMessage('Background subtraction, %s_flt.fits:\n\n  %s' %(exp, coeff_str))
        
    return models
