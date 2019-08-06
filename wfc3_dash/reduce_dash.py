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

    Lacosmicx (used to identify cosmic rays): https://github.com/cmccully/lacosmicx

References
----------

    - Ivelina G. Momcheva, Pieter G. van Dokkum, Arjen van der Wel,
    Gabriel B. Brammer, John MacKenty, Erica J. Nelson, Joel Leja,
    Adam Muzzin, and Marijn Franx. 2017 January. 
    doi:10.1088/1538-3873/129/971/015004

"""

from astropy.io import fits
import numpy as np
from drizzlepac import tweakreg
from drizzlepac import astrodrizzle
import stwcs
import os
from utils import get_flat
from utils import get_IDCtable
from glob import glob #Needed for gaia alignment
from stsci.tools import teal #Needed for gaia alignment
from stwcs import updatewcs #Needed for gaia alignment
from astropy.convolution import Gaussian2DKernel #Needed for create_seg_map
from astropy.stats import gaussian_fwhm_to_sigma #Needed for create_seg_map
from astropy.io import ascii
from astropy.table import Table, Column, MaskedColumn
from photutils import detect_sources #Needed for create_seg_map
from photutils import detect_threshold #Needed for create_seg_map
from photutils import source_properties
import lacosmicx #for fix_cosmic_rays
import os

class DashData(object):
    
    def __init__(self,file_name=None, flt_file_name=None):
        '''
		The init method performs a series of tests to make sure that the file fed to the DashData class is a valid IMA file with units in e/s. Will also add root attribute to DashData class.

		Paramaters
		---------------
		self : object
			DashData object created from an individual IMA file
		file_name : fits file
			Name of IMA file that is fed to DashData class
        flt_file_name : fits file
            Name of FLT file that is fed to DashData class

		Outputs
		-----------
		N/A
		'''

        if '_ima.fits' not in file_name: 
            raise Exception('Input needs to be an IMA file.')
        else:
            self.file_name = file_name
        #First test whether the file exists
            try:
                self.ima_file = fits.open(self.file_name)
            ### If it does - test its content

            ###### Test whether it is a wfc3/ir image
                if ( (self.ima_file[0].header['INSTRUME'].strip() == 'WFC3') & (self.ima_file[0].header['DETECTOR'].strip() == 'IR') ) == False:
                    raise Exception('This observation was not performed with WFC3/IR, instrument is set to: {}, and detector is set to: {}'.format(
                                    self.ima_file[0].header['INSTRUME'],self.ima_file[0].header['DETECTOR'] ) )

                ###### Test whether it has more than one science extension (FLTs have only 1)
                nsci = [ext.header['EXTNAME '] for ext in self.ima_file[1:]].count('SCI')
                if nsci == 1:
                    raise Exception('This file has only one science extension, it cannot be a WFC3/IR ima')

                ###### Test that this file is NOT a RAW
                calib_keys = ['DQICORR','ZSIGCORR','ZOFFCORR','DARKCORR','BLEVCORR','NLINCORR','FLATCORR','CRCORR','UNITCORR','PHOTCORR','RPTCORR','DRIZCORR']
                performed = 0
                for ck in calib_keys:
                    if self.ima_file[0].header[ck] == 'COMPLETE':
                        performed = performed + 1
                if performed == 0:
                    raise Exception('This file looks like a RAW file')     

                ###### Test that the units of the individual images are e/s
                
                bu = self.ima_file[1].header['BUNIT']
                if  bu != 'ELECTRONS/S':
                    uc = self.ima_file[0].header['UNITCORR']
                    fc = self.ima_file[0].header['FLATCORR']
                    raise Exception('BUNIT in the "SCI" extensions of this file is set to "{}", but should be set to "ELECTRONS/S"\n'
                                    'This is a consequence of UNITCORR set to "{}" and FLATCORR set to "{}".\n'
                                    'Please rerun calwf3 on this file after setting both UNITCORR and FLATCORR to "PERFORM" in the 0-th extension header'.format(bu,uc,fc))
                    
                

            except IOError:
                print('Cannot read file.')
        
        self.flt_file_name = flt_file_name  
        self.root = self.file_name.split('/')[-1].split('_ima')[0]
    
    def align(self, subtract_background = True, 
              align_method = None, ref_catalog = None, 
              create_diff_source_lists=True,
              updatehdr=True, updateWCS=True, wcsname = 'DASH', 
              threshold = 50., cw = 3.5, 
              searchrad=20., astrodriz=True, 
              cat_file='diff_catfile.cat',
              drz_output=None,):

        '''
        Aligns new FLT's to reference catalog.

        Parameters
        ----------
        self : object
            DashData object created from an individual IMA file.
        subtract_background : bool, optional
            If True, runs subtract_background_reads functions.
        align_method : string, optional
            Defines alignment method to be used. Default is None (input files will align to each other).
        ref_catalog : cat file, optional
            Defines reference image that will be referenced for CATALOG alignment method.
        create_diff_source_lists : bool, optional
            Specifies whether or not to create a segmentation image and source list from the difference files
        updatehdr : bool, optional
            Specifies whether to update the headers after aligning during TweakReg. Default is True.
        updateWCS : bool, optional
            Specifies whether to update the WCS after aligning during TweakReg. Default is True.
        wcsname : str, optional (as long as name you choose doesn't already exist)
            Specifies name of WCS. Default is 'DASH'
        threshold : float, optional
        cw : float, optional
        searchrad : float, optional
            Radius (in pixels) that TweakReg will search around sources to find matches. Default is 20 pixels.
        astrodriz : bool, optional
            Specifies whether to drizzle images together using Astrodrizzle. Default is True
        cat_file : str, optional
            Name of catfile to be used to align sources in TweakReg. Default is the catfile created by setting create_diff_source_lists to True, diff_catfile.cat
        drz_output : str, optional
            Name of output file after drizzling using AstroDrizzle. Default is the root name of the original IMA.


        Outputs
        -------
        Shifts file : txt file
            File containing shifts during TweakReg
        WCS Shifts file : fits file
            File containing WCS shifts during TweakReg
        Residual plot : png
            Plot showing residuals from TweakReg alignment
        Vector plot : png
            Plot showing vector poitnings from TweakReg alignment
        Aligned FLT's : fits file
            Same diff files created in split_ima that have been aligned using TweakReg
        Drizzled Image : fits file
            Setting astrodriz to True will output a drizzled image form a single exposure (produced only if astrodriz is set to True)
        '''

        #Set name for output drizzled image to the rootname of the original IMA if it is not specified
        if drz_output is None:
            drz_output=self.root

        input_images = sorted(glob('diff/{}_*_diff.fits'.format(self.root)))

        outshifts = 'shifts_{}.txt'.format(self.root)
        outwcs = 'shifts_{}_wcs.fits'.format(self.root)
        
        if subtract_background:
            self.subtract_background_reads()
        
        #Align images to a catalog
        if align_method == 'CATALOG': 
            if (ref_catalog is not None):

                ##Create source list and segmentation maps based on difference files
                if create_diff_source_lists is True:
                    diffpath = os.path.dirname(os.path.abspath('diff/{}_*_diff.fits'.format(self.root)))
                    cat_images=sorted([os.path.basename(x) for x in glob('diff/{}_*_diff.fits'.format(self.root))])

                    sc_diff_files = [diffpath + '/' + s for s in cat_images]

                    self.diff_seg_map(cat_images=sc_diff_files)

                wcs = list(map(updatewcs.updatewcs, input_images))

                teal.unlearn('tweakreg')
                teal.unlearn('imagefindpars')

                tweakreg.TweakReg(input_images, 
                                  refcat=ref_catalog,
                                  catfile=cat_file,
                                  xcol=2,
                                  ycol=3,
                                  updatehdr=updatehdr, 
                                  updatewcs=updateWCS,
                                  wcsname=wcsname,
                                  verbose=True,
                                  imagefindcfg={'threshold': threshold, 'conv_width': cw},
                                  searchrad=searchrad,
                                  searchunits = 'pixels',
                                  shiftfile=True, 
                                  outshifts=outshifts,
                                  outwcs=outwcs,
                                  interactive=False,
                                  wcsname='DASH',  
                                  fitgeometry='rscale',
                                  minobj=5)

            else:
                
                raise Exception('Need to specify reference catalog and reference image.')
            
        
        #Align images to the first image        
        else:
            
                wcs = list(map(updatewcs.updatewcs, input_images))

                teal.unlearn('tweakreg')
                teal.unlearn('imagefindpars')

                tweakreg.TweakReg(input_images, 
                                  updatehdr=updatehdr, 
                                  updatewcs=updateWCS,
                                  writecat=True,
                                  verbose=True,
                                  imagefindcfg={'threshold':threshold,'conv_width':cw},
                                  refcat=ref_catalog,
                                  clean=False, 
                                  interactive=False,
                                  see2dplot=False,
                                  shiftfile=True, 
                                  outshifts=outshifts,
                                  outwcs=outwcs,
                                  wcsname=wcsname, 
                                  headerlet = False,
                                  minobj = 5, 
                                  searchrad = searchrad, 
                                  searchunits = 'pixels', 
                                  use2dhist = True,
                                  reusename=True,
                                  fitgeometry='rscale',
                                  nclip=3,
                                  sigma=3.0,
                                  clobber=True,
                                  dqbits=0)
            
                pass

        #Drizzle the images together
        if astrodriz is True:    

            #Do not have drizzle take 256 flags into account
            no_tfs = 2,4,8,16,32,64,128,512,2048,4096,8192,16384
    
            astrodrizzle.AstroDrizzle(input_images, 
                output=drz_output,
                clean=False, 
                final_pixfrac=1.0, 
                context=False,
                resetbits=0, 
                preserve=False, 
                driz_cr_snr='8.0 5.0', 
                driz_cr_scale = '2.5 0.7',
                driz_sep_bits=no_tfs,
                final_bits=no_tfs)    

    def create_seg_map(self):
        '''
        Creates segmentation map, from original FLT file, that is used in background subtraction and to fix cosmic rays.

        Parameters
        ----------
        self : object
            DashData object created from an individual IMA file.

        Output
        ------
        Segmentation Image : fits file
            Segmentation map
        Source List : .dat file
            List of sources and their properties
        '''
    
        flt = fits.open(self.flt_file_name)
        data = flt[1].data

        threshold = detect_threshold(data, snr=3.)

        sigma = 3.0 * gaussian_fwhm_to_sigma    # FWHM = 3.
        kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
        kernel.normalize()
        segm = detect_sources(data, threshold, npixels=10, filter_kernel=kernel)

        hdu = fits.PrimaryHDU(segm.data)
        hdu.writeto(('{}_seg.fits').format(self.root), clobber=True)

        # Create source list
        cat = source_properties(data, segm)

        tbl = cat.to_table()
        tbl['xcentroid'].info.format = '.2f'
        tbl['ycentroid'].info.format = '.2f'
        tbl['cxx'].info.format = '.2f'
        tbl['cxy'].info.format = '.2f'
        tbl['cyy'].info.format = '.2f'

        ascii.write(tbl, '{}_source_list.dat'.format(self.root))

    def diff_seg_map(self, cat_images=None, remove_column_names=True, snr=1.0, sig=6.0, npixels=5):
        '''
        Creates segmentation image and source list from difference files.

        Parameters
        ----------
        self : object
            DashData object created from an individual IMA file
        cat_images : list, str
            List of difference files with full path name
        remove_column_names : bool
            Specifies whether to remove the header from the source lists so TweakReg can read them
        snr : float
        sig : float
        npixels : float

        Outputs
        -------
        Segmentation Image : fits file
            Segmentation map
        Source List : .dat file
            List of sources and their properties
        '''

        input_images = sorted(glob('diff/{}_*_diff.fits'.format(self.root)))

        for index, diff in enumerate(input_images, start=1):

            diff=fits.open(diff)
            data = diff[1].data

            threshold = detect_threshold(data, snr=snr)

            sigma = sig * gaussian_fwhm_to_sigma 
            kernel = Gaussian2DKernel(sigma, x_size=sig, y_size=sig)
            kernel.normalize()
            segm = detect_sources(data, threshold, npixels=npixels, filter_kernel=kernel)

            hdu = fits.PrimaryHDU(segm.data)
            hdu.writeto(('{}_{:02d}_diff_seg.fits').format(self.root, index), clobber=True)

            # Create source list
            cat = source_properties(data, segm)

            tbl = cat.to_table()
            tbl['xcentroid'].info.format = '.2f'
            tbl['ycentroid'].info.format = '.2f'
            tbl['cxx'].info.format = '.2f'
            tbl['cxy'].info.format = '.2f'
            tbl['cyy'].info.format = '.2f'

            ascii.write(tbl, '{}_{:02d}_diff_source_list.dat'.format(self.root, index))

            if remove_column_names is True:

                #Remove headers to source lists so tweakreg can read them
                n = 1
                nfirstlines = []

                with open('{}_{:02d}_diff_source_list.dat'.format(self.root, index)) as f, open("temp_sl.dat", "w") as out:
                    for x in range(n):
                        nfirstlines.append(next(f))
                    for line in f:
                        out.write(line)

                os.remove('{}_{:02d}_diff_source_list.dat'.format(self.root, index))
                os.rename("temp_sl.dat", '{}_{:02d}_diff_source_list.dat'.format(self.root, index))
        
        if cat_images is not None:
            x=np.array(cat_images)
            y=np.array(sorted(glob(('{}_*_diff_source_list.dat').format(self.root))))
            catdata = Table([x, y], names=['Diff File', 'Source List'])
            ascii.write(catdata, 'diff_catfile.cat')
        else:
            raise Exception('Need to input list of difference files in order to make source list. List should include full path.')


    def fix_cosmic_rays(self, rm_custom=False, flag=None):
        '''
        Resets cosmic rays within the seg maps of objects and uses L.A.Cosmic to find them again.

        Parameters
        ----------
        self : object
            DashData object created from an individual IMA file.
        rm_custom : bool
            Specifies whether or not the user would like to remove custom flags within the boundaries of sources, as defined by the segmentation map created from the original FLT.
        flag : int
            Specifies flag the user would like the remove within the boundaries of sources.


        Output
        ------
        Fixed for cosmic rays diff files : fits
            Same diff files created in split_ima that have now been corrected for cosmic ray errors.
        '''

        asn_exposures = sorted(glob('diff/' + self.root + '_*_diff.fits'))

        seg = fits.open('{}_seg.fits'.format(self.root))
        seg_data = np.cast[np.float32](seg[0].data)

        flt_full = fits.open(self.flt_file_name)
        flt_full_wcs = stwcs.wcsutil.HSTWCS(flt_full, ext=1)

        EXPTIME = flt_full[0].header['EXPTIME']

        crmask, clean = lacosmicx.lacosmicx(flt_full[1].data, gain=1.0, readnoise=20., 
                                    objlim = 15.0, 
                                    pssl = 0., 
                                    verbose=True)

        yi, xi = np.indices((1014,1014))

        for exp in asn_exposures:

            flt = fits.open(exp, mode = 'update')

            flagged_stars = ((flt['DQ'].data & 4096) > 0) & (seg_data > 0)
            flt['DQ'].data[flagged_stars] -= 4096

            new_cr = (crmask == 1) & ((flt['DQ'].data & 4096) == 0) & ((seg_data == 0) | ((seg_data > 0) & (flt['SCI'].data < 1.)))  & (xi > 915) & (yi < 295)

            flt['DQ'].data[new_cr] += 4096

            flt.flush()

        if rm_custom is True:

            if flag is not None:

                for exp in asn_exposures:

                    flt = fits.open(exp, mode = 'update')

                    flagged_stars = ((flt['DQ'].data & flag) > 0) & (seg_data > 0)
                    flt['DQ'].data[flagged_stars] -= flag

                    new_cr = (crmask == 1) & ((flt['DQ'].data & flag) == 0) & ((seg_data == 0) | ((seg_data > 0) & (flt['SCI'].data < 1.)))  & (xi > 915) & (yi < 295)

                    flt['DQ'].data[new_cr] += flag

                    flt.flush()

            else:

                raise Exception('Must specify which flags to remove.')

    def make_pointing_asn(self):
        """ 
        Makes a new association table for the reads extracted from a given IMA.
        Parameters
        ----------
        self : object
            DashData object created from an individual IMA file. 
        Outputs
        ----------
        ASN files : fits file
            Fits file formatted as an association file that holds 
            the names of the difference files created by split_ima
            and the root name of the individual IMA file.
        """
        
        asn_filename = 'diff/{}_asn.fits'.format(self.root)
        asn_list = self.diff_files_list.copy()
        asn_list.append(self.root)

        # Create Primary HDU:
        hdr = fits.Header()
        hdr['FILENAME'] = asn_filename
        hdr['FILETYPE'] = 'ASN_TABLE'
        hdr['ASN_ID'] = self.root
        hdr['ASN_TABLE'] = asn_filename
        hdr['COMMENT'] = "This association table is for the read differences for the IMA."
        primary_hdu = fits.PrimaryHDU(header=hdr)

        # Create the information in the asn file
        num_mem = len(asn_list)

        asn_mem_names = np.array(asn_list)
        asn_mem_types =  np.full(num_mem,'EXP-DTH',dtype=np.chararray)
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
             
    def run_reduction(self):
        ''' 
        wrapper
        '''

        pass              

    def split_ima(self):
        ''' 
        Will create individual files for the difference between 
        adjacent reads of a IMA file. Will also add more attributes 
        to the DashData object. 
        Parameters
        ----------
        self : object
            DashData object created from an individual IMA file. 
        Outputs
        ----------
        N files : fits files
            Fits files of the difference between adjacent IMA reads.
        
        '''
        FLAT = fits.open(get_flat(self.file_name))
        IDCtable = fits.open(get_IDCtable(self.file_name))
        
        NSAMP = self.ima_file[0].header['NSAMP']
        shape = self.ima_file['SCI',1].shape
    
        cube = np.zeros((NSAMP, shape[0], shape[1]), dtype='float32')
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
        
        self.readnoise_2D = np.zeros((1024,1024), dtype='float32')
        self.readnoise_2D[512: ,0:512] += self.ima_file[0].header['READNSEA']
        self.readnoise_2D[0:512,0:512] += self.ima_file[0].header['READNSEB']
        self.readnoise_2D[0:512, 512:] += self.ima_file[0].header['READNSEC']
        self.readnoise_2D[512: , 512:] += self.ima_file[0].header['READNSED']
        self.readnoise_2D = self.readnoise_2D**2
        
        self.diff_files_list = []
        for j in range(1, NSAMP-1):

            hdu0 = fits.PrimaryHDU(header=self.ima_file[0].header)
            hdu0.header['EXPTIME'] = dt[j]
            hdu0.header['IMA2FLT'] = (1, 'FLT {} extracted from IMA file'.format(j)) 
            hdu0.header['NEXTEND'] = 5
            hdu0.header['OBSMODE'] = 'ACCUM'
            hdu0.header['NSAMP'] = 1

            
            # NOT SURE THE HEADERS I AM GIVING IT HERE ARE OK
            science_data = diff[j,:,:]/dt[j]
            hdu1 = fits.ImageHDU(data = science_data[5:-5,5:-5], header = self.ima_file['SCI',NSAMP-j-1].header, name='SCI')
            hdu1.header['EXTVER'] = 1
            hdu1.header['ROOTNAME'] = '{}_{:02d}'.format(self.root,j)

            var = 2*self.readnoise_2D + science_data*FLAT['SCI'].data*dt[j]
            err = np.sqrt(var)/dt[j]
                    
            hdu2 = fits.ImageHDU(data = err[5:-5,5:-5], header = self.ima_file['ERR',NSAMP-j-1].header, name='ERR')
            hdu2.header['EXTVER'] = 1
            
            hdu3 = fits.ImageHDU(data = dq[j+1][5:-5,5:-5], header = self.ima_file['DQ',NSAMP-j-1].header, name='DQ')
            hdu3.header['EXTVER'] = 1
            
            #hdu3.data[BP_MASK == 1] += 4
            # trun the 8192 cosmic ray flag to the standard 3096
            hdu3.data[(hdu3.data & 8192) > 0] -= 4096
            # remove the 32 flag, these are not consistently bad
            hdu3.data[(hdu3.data & 32) > 0] -= 32
            
            hdu4 = fits.ImageHDU(data = (np.zeros((1014,1014), dtype=np.int16) + 1), header = self.ima_file['SAMP',NSAMP-j-1].header, name = 'SAMP')
            hdu5 = fits.ImageHDU(np.zeros((1014,1014)) + dt[j], header = self.ima_file['TIME',NSAMP-j-1].header, name = 'TIME')
            hdu4.header['EXTVER'] = 1
            hdu5.header['EXTVER'] = 1

            hdu = fits.HDUList([hdu0,hdu1,hdu2,hdu3,hdu4,hdu5])
            print('Writing {}_{:02d}_diff.fits'.format(self.root,j))

            self.hdu = hdu

            if not os.path.exists('diff'):
                os.mkdir('diff')

            
            hdu.writeto('diff/{}_{:02d}_diff.fits'.format(self.root,j), overwrite=True)
            
            self.diff_files_list.append('diff/{}_{:02d}'.format(self.root,j))
        
    def subtract_background_reads(self, subtract=True, reset_stars_dq=False):
        '''
        Performs median background subtraction for each individual difference file.
        Uses the DRZ and SEG images produced in FLT background subtraction.
        Parameters
        ----------
        self : object
            DashData object created from an individual IMA file.
        subtract : bool
            Does not subtract the background by default, but still writes it to the header.
            Set to True to subtract background.
        reset_stars_dq : bool
            Set to True to reset cosmic rays within objects to 0 because the centers of stars are flagged.
        Outputs
        -------
        Background Subtracted N files : fits files
            Fits files of the difference between adjacent IMA reads that have been background subtracted.
        '''

        ### I think this should subtract the background on only one FLT
        ### But currently loops over all of them

        seg = fits.open('{}_seg.fits'.format(self.root))    
        seg_data = np.cast[np.float32](seg[0].data)
        
        yi, xi = np.indices((1014,1014))
        
        
        self.bg_models = []
        
        for ii, exp in enumerate(self.diff_files_list):
        
            diff = fits.open('{}_diff.fits'.format(exp), mode='update')
            diff_wcs = stwcs.wcsutil.HSTWCS(diff, ext=1)
                
            mask = (seg_data == 0) & (diff['DQ'].data == 0) & (diff[1].data > -1) & (xi > 10) & (yi > 10) & (xi < 1004) & (yi < 1004)
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
                    
                    
def main(ima_file_name = None, flt_file_name = None, 
         align_method = None, ref_catalog = None, 
         drz_output=None, subtract_background = False, 
         wcsname = 'DASH', threshold = 50., cw = 3.5, 
         updatehdr=True, updateWCS=True, 
         searchrad=20., 
         astrodriz=True):

    '''
    Runs entire DashData pipeline under a single function.

    Parameters
    ----------
    ima_file_name : str
        File name of ima file.
    flt_file_name : str
        File name of FLT file.
    align_method : str, optional
        Method to align difference files using TweakReg. Default is None, which aligns reads to the first read.
        Setting align_method equal to 'CATALOG' will align the reads to a catalog.
    ref_catalog : str, optional
        Catalog to be aligned to if using CATALOG align method.
    drz_output : str, optional
        Name to call final drizzled output image ( + '_drz_sci.fits' ). Default is root_name + '_drz_sci.fits'.
    subtract_background : bool
        Determines whether background is subtracted or not during align function. Default is False since during this main function,
        the background is subtracted separately.
    wcsname : str
        Name for WCS during TweakReg.
    threshold : float
        qwerty
    cw : float
        qwerty
    updatehdr : bool
        Determines whether to update headers during TweakReg. Default is True.
    updateWCS : bool
        Determines whether to update WCS information during TweakReg. Default is True.
    searchrad : float
        qwerty
    astrodriz : bool
        Determines whether or not to run astrodrizzle. Default is True.

    Outputs
    -------
    N files : fits files
        Fits files of the difference between adjacent IMA reads.
    N files : fits files
        Fits files of the difference between adjacent IMA reads.
    Shifts file : txt file
        File containing shifts during TweakReg
    WCS Shifts file : fits filr
        File containing WCS shifts during TweakReg
    Drizzled science image : fits
        Drizzled science image from one exposure reduced using the DASH pipeline.
    Weighted science image : fits
        Weighted drizzled science image from one exposure reduced using the DASH pipeline.
    '''
    
    myDash = DashData(ima_file_name, flt_file_name)
    myDash.split_ima()
    myDash.create_seg_map()
    myDash.subtract_background_reads()
    myDash.fix_cosmic_rays()
    myDash.align(align_method = align_method, ref_catalog = ref_catalog, drz_output=drz_output, 
                 subtract_background = subtract_background, 
                 wcsname = wcsname, threshold = threshold, cw = cw, 
                 updatehdr=updatehdr, updateWCS=updateWCS, 
                 searchrad=searchrad, 
                 astrodriz=astrodriz)


