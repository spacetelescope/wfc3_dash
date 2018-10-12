#! /usr/bin/env python

""" Reduces DASH/IR data. 

This script serves as a pipeline to reduce DASH/IR data. The 
FLT/IMA files available from MAST are inappropiate to use to 
perform science upon (Momcheva et al. 2017). This script will 
create an aligned image that is the combination of the 
differences between adjacent reads. This file is appropriate 
to use for science work. 


The following output products are created by this script from 
each IMA file:

1. Individual fits files of differences between IMA reads with 
    thier own headers. 
2. A master ASN file of those read fits files for a single IMA. 
3. An aligned version of the IMA read fits files. 

Authors
-------
    Iva Momcheva 2018
    Catherine Martlin 2018
    Mario Gennaro 2018

Use
---
    This script is intended to be executed via the command line 
    as such:
    ::

        python reduce_dash.py [-f|--file]

    ``-f --file`` - The IMA file name/path. 

Notes
-----

References
----------

    - Ivelina G. Momcheva, Pieter G. van Dokkum, Arjen van der Wel,
    Gabriel B. Brammer, John MacKenty, Erica J. Nelson, Joel Leja,
    Adam Muzzin, and Marijn Franx. 2017 January. 
    doi:10.1088/1538-3873/129/971/015004

"""
import numpy as np

from astropy.io import fits
import numpy as np
from urllib.request import urlretrieve 
from drizzlepac import tweakreg
import stwcs
import os

PATH = 'data_download/mastDownload/HST/'


class DashData(object):
    
    def __init__(self,file_name=None):
        
        if '_ima.fits' not in file_name: 
            raise Exception('Input needs to be an IMA file.')
        else:
            self.file_name = file_name
	    #First test whether the file exists
            try:
                self.ima_file = fits.open(self.file_name)
		### If it does test its content

		###### Test whether it is a wfc3/ir image
		if ~( (INSTRUME == 'WFC3  ') & (DETECTOR == 'IR  ') ):
		    raise Exception('This observation was not performed with WFC3/IR')

		###### Test whether it has more than one science extension (FLTs have only 1)
	        nsci = [ext.header['EXTNAME '] for ext in self.ima_file[1:]].count('SCI')
		if nsci == 1:
		    raise Exception('This file has only one science extension, it cannot be a WFC3/IR ima')


		calib_keys = ['DQICORR','ZSIGCORR','ZOFFCORR','DARKCORR','BLEVCORR','NLINCORR','FLATCORR','CRCORR','UNITCORR','PHOTCORR','RPTCORR','DRIZCORR']
		performed = 0
		for ck in calib_keys:
		    if self.ima_file[ck] == 'COMPLETE':
			performed = performed + 1
		if performed == 0:
		    raise Exception('This file looks like a RAW file')     
		
            except IOError:
                print('Cannot read file.')
                
        self.root = self.file_name.split('/')[-1].split('_ima')[0]

        ### Need to check header and make sure that this is actually a gyro exposure

	
        

    def split_ima(self):
        
        FLAT = fits.open(get_flat(self.file_name))
        
        NSAMP = self.ima_file[0].header['NSAMP']
        shape = self.ima_file['SCI',1].shape
    
        cube = np.zeros((NSAMP, shape[0], shape[1]))
        dq = np.zeros((NSAMP, shape[0], shape[1]), dtype=np.int)
        time = np.zeros(NSAMP)
        
        for i in range(NSAMP):
            cube[NSAMP-1-i, :, :] = self.ima_file['SCI',i+1].data*self.ima_file['TIME',i+1].header['PIXVALUE']
            dq[NSAMP-1-i, :, :] = self.ima_file['DQ',i+1].data
            time[NSAMP-1-i] = self.ima_file['TIME',i+1].header['PIXVALUE']
        
        diff = np.diff(cube, axis=0)
        self.diff = diff[1:]
        dt = np.diff(time)
        self.dt = dt[1:]
        
        self.dq = dq[1:]
        
        self.readnoise_2D = np.zeros((1024,1024))
        self.readnoise_2D[512: ,0:512] += self.ima_file[0].header['READNSEA']
        self.readnoise_2D[0:512,0:512] += self.ima_file[0].header['READNSEB']
        self.readnoise_2D[0:512, 512:] += self.ima_file[0].header['READNSEC']
        self.readnoise_2D[512: , 512:] += self.ima_file[0].header['READNSED']
        self.readnoise_2D = self.readnoise_2D**2
        
        self.file_list = []
        for j in range(1, NSAMP-1):

            hdu0 = fits.PrimaryHDU(header=self.ima_file[0].header)
            hdu0.header['EXPTIME'] = dt[j]
            hdu0.header['IMA2FLT'] = (1, 'FLT {} extracted from IMA file'.format(j)) 
            
            # NOT SURE THE HEADERS I AM GIVING IT HERE ARE OK
            science_data = diff[j,:,:]/dt[j]
            hdu1 = fits.ImageHDU(data = science_data[5:-5,5:-5], header = self.ima_file['SCI',NSAMP-j-1].header, name='SCI')

            var = 2*self.readnoise_2D + science_data*FLAT['SCI'].data
            err = np.sqrt(var)/dt[j]
                    
            hdu2 = fits.ImageHDU(data = err[5:-5,5:-5], header = self.ima_file['ERR',NSAMP-j-1].header, name='ERR')
            
            hdu3 = fits.ImageHDU(data = dq[j+1][5:-5,5:-5], header = self.ima_file['DQ',NSAMP-j-1].header, name='DQ')
            
            #hdu3.data[BP_MASK == 1] += 4
            # trun the 8192 cosmic ray flag to the standard 3096
            hdu3.data[(hdu3.data & 8192) > 0] -= 4096
            # remove the 32 flag, these are not consistently bad
            hdu3.data[(hdu3.data & 32) > 0] -= 32
            
            hdu4 = fits.ImageHDU(data = (np.zeros((1014,1014)) + 1.), name = 'SAMP')
            hdu5 = fits.ImageHDU(np.zeros((1014,1014)) + dt[j], name = 'TIME')
            
            hdu = fits.HDUList([hdu0,hdu1,hdu2,hdu3,hdu4,hdu5])
            print('Writing {}_{:02d}_diff.fits'.format(self.root,j))
            hdu.writeto('{}_{:02d}_diff.fits'.format(self.root,j), overwrite=True)
            
            self.file_list.append('{}_{:02d}'.format(self.root,j))
        
    def subtract_background_reads(self, subtract=True, reset_stars_dq=False):
        
        ### I think this should subtract the background on only one FLT
        ### But currently loops over all of them
        
        # read in DRZ and SEG images produced for the FLT
        # should these be saved in the class variables?
        
        ### UNCOMMENT THESE ONCE AVAILABLE
        #ref = fits.open('{}_drz_sci.fits'.format(self.root))
        #ref_wcs = stwcs.wcsutil.HSTWCS(ref, ext=0)

        #seg = fits.open('{}_drz_seg.fits'.format(self.root))    
        #seg_data = np.cast[np.float32](seg[0].data)
        
        yi, xi = np.indices((1014,1014))
        
        
        self.bg_models = []
        
        for ii, exp in enumerate(self.file_list):
        
            diff = fits.open('{}_diff.fits'.format(exp), mode='update')
            diff_wcs = stwcs.wcsutil.HSTWCS(diff, ext=1)
            
            if ii == 0:
                ### Only done for the first one because it is just used as a mask
                print('Segmentation image: {}_blot.fits'.format(exp))
                
                ### UNCOMMENT THESE ONCE THE SEG MAP IS AVAILABLE
                #blotted_seg = astrodrizzle.ablot.do_blot(seg_data, ref_wcs, diff_wcs, 1, coeffs=True, interp='nearest', sinscl=1.0, stepsize=10, wcsmap=None)         
                blotted_seg = np.zeros((1014,1014))
                
            mask = (blotted_seg == 0) & (diff['DQ'].data == 0) & (diff[1].data > -1) & (xi > 10) & (yi > 10) & (xi < 1004) & (yi < 1004)
            mask &= (diff[1].data < 5*np.median(diff[1].data[mask]))
            data_range = np.percentile(diff[1].data[mask], [2.5, 97.5])
            mask &= (diff[1].data >= data_range[0]) & (diff[1].data <= data_range[1])
            data_range = np.percentile(diff[2].data[mask], [0.05, 99.5])
            mask &= (diff[2].data >= data_range[0]) & (diff[2].data <= data_range[1])
        
            sky_level = np.median(diff[1].data[mask])
            model = diff[1].data*0. + sky_level
        
            # add header keywords of the fit components

            diff[1].header['MDRIZSKY'] =  sky_level 
            if not subtract:
                if 'BG_SUB' not in (key for key in diff[1].header.keys()):
                    flt[1].header['BG_SUB'] =  'No'
            else:
                if 'BG_SUB' not in (key for key in diff[1].header.keys()):
                    diff[1].data -= model           
                    diff[1].header['BG_SUB'] =  'Yes'
                else:
                    print('Keyword BG_SUB already set to {}. Skipping background subtraction.'.format(diff[1].header['BG_SUB']))
                
        
            if reset_stars_dq:
                flagged_stars = ((diff['DQ'].data & 4096) > 0) & (blotted_seg > 0)
                diff['DQ'].data[flagged_stars] -= 4096
        
            diff.flush()
            print('Background subtraction, {}_diff.fits:  {}'.format(exp, sky_level))
            
            ### Should these background subtracted reads substitute the ones saved in split_ima?
                    
    def subtract_background_flt():
        
        pass
        
    def fix_cosmic_rays():
        
        pass
        
    def make_read_catalog():
        
        pass
        
    def align_read():
        
        pass
        
    def align(self, align_method = 'CATALOG', ref_catalog = None, ref_image = None, subtract_background = True):
        
        outshifts = 'shifts_{}.txt'.format(self.root)
        outwcs = 'shifts_{}_wcs.fits'.format(self.root)
        
        if subtract_background:
            self.subtract_background_reads()
        
        ### What other alignment methods are there? TWEAKREG, GAIA?
        
        if align_method == 'CATALOG': 
            if (ref_catalog is not None) and (ref_image is not None):
                
                drizzlepac.tweakreg.TweakReg(self.root+'_flt.fits', 
                    refimage=ref_image, # reference image
                    refcat = ref_catalog, # reference catalog
                    updatehdr = True, 
                    updatewcs = True, 
                    writecat = False, 
                    clean = True, 
                    verbose = True, 
                    runfile = 'tweakreg.log', 
                    wcsname = 'TWEAK', 
                    headerlet = False, 
                    shiftfile = True, 
                    outshifts = outshifts, 
                    outwcs = outwcs, 
                    refxcol = 5, 
                    refycol = 6, 
                    refxyunits = 'degrees', 
                    minobj = 5, 
                    searchrad = 1000.0, 
                    searchunits = 'pixels', 
                    use2dhist = True, 
                    see2dplot = False, 
                    separation = 0.5, 
                    tolerance = 1.0, 
                    xoffset = 0.0, 
                    yoffset = 0.0, 
                    fitgeometry = 'shift', 
                    interactive=False, 
                    nclip = 3, 
                    sigma = 3.0, 
                    clobber=True) 
            else:
                
                ### what do you do if there's no catalog - throw an error?
                
                pass
                
        else:
            ### Needs to actually do tweakreg?
            
            pass
            
    
        drizzlepac.astrodrizzle.AstroDrizzle(root+'_flt.fits', 
            clean=False, 
            final_pixfrac=1.0, 
            context=False, 
            final_bits=576, 
            resetbits=0, 
            preserve=False, 
            driz_cr_snr='8.0 5.0', 
            driz_cr_scale = '2.5 0.7', 
            wcskey= 'TWEAK')
              
    def coadd_reads():
        
        pass

    def make_pointing_asn(self):
        """ Makes a new association table for the reads extracted from a given IMA.

        """
        asn_filename = '{}_asn.fits'.format(self.root)
        file_list = self.file_list
        asn_list = file_list.append(self.root)

        # Create Primary HDU:
        hdr = fits.Header()
        hdr['FILENAME'] = "'" + asn_filename + "'"
        hdr['FILETYPE'] = 'ASN_TABLE'
        hdr['ASN_ID'] = "'" + self.root "'"
        hdr['ASN_TABLE'] = "'" + asn_filename + "'"
        hdr['COMMENT'] = "This association table is for the read differences for the IMA."
        primary_hdu = fits.PrimaryHDU(header=hdr)

        # Create the information in the asn file
        num_mem = len(asn_list)

        asn_mem_names = np.array(asn_list)
        asn_mem_types =  (np.full(num_mem,'EXP-DTH'))
        asn_mem_types[-1] = 'PROD-DTH'
        asn_mem_prsnt = np.ones(num_mem, dtype=np.bool_)
        asn_mem_prsnt[-1] = 0

        hdu_data = fits.BinTableHDU().from_columns([fits.Column(name='MEMNAME', format='14A', array=asn_mem_names), 
                    fits.Column(name='MEMTYPE', format='14A', array=asn_mem_types), 
                    fits.Column(name='MEMPRSNT', format='L', array=asn_mem_prsnt)])

        # Create the final asn file
        hdu = fits.HDUList([primary_hdu, hdu_data])

        if 'EXTEND' not in hdu[0].header.keys():
            hdu[0].header.update('EXTEND', True, after='NAXIS')
        
        hdu.writeto(asn_filename, overwrite=True)

        # Create property of the object that is the asn filename.
        self.asn_filename = asn_filename 

    def run_reduction():
        pass

def main():
    ''' 
    Main function of reduce_dash. 
    
    Parameters
    ----------
     : 
        
    Outputs
    ----------
     : 
         
    '''
    
    pass


