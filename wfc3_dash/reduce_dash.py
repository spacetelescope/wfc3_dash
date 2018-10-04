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
from urllib.request import urlretrieve 

BP_MASK = fits.open('')[0].data
PATH = 'data_download/mastDownload/HST/'


class DashData:
    
    def __init__(self,file):
        
        if '_ima.fits' not in file:
            raise Exception('Input needs to be an IMA file.')
        else:
            self.file = file
            try:
                ima_file = fits.open(self.file)
            except IOError:
                print('Cannot read file.')

        ### Need to check header and make sure that this is actually a gyro exposure
        
    
    def split_ima():

        pass
        
    def subtract_background_reads():
        
        pass
        
    def fix_cosmic_rays():
        
        pass
        
    def align_reads():
        
        pass
        
    def coadd_reads():
        
        pass

    def make_pointing_asn(self):
        """ Makes a new association table for the reads extracted from a given IMA.

        """
        root = self.root
        asn_filename = root +'_asn.fits'



        # Create Primary HDU:
        hdr = fits.Header()
        hdr['FILENAME'] = "'" + asn_filename + "'"
        hdr['FILETYPE'] = 'ASN_TABLE'
        hdr['ASN_ID'] = "'" + root "'"
        hdr['ASN_TABLE'] = "'" + asn_filename + "'"
        hdr['COMMENT'] = "This association table is for the read differences for the IMA."
        primary_hdu = fits.PrimaryHDU(header=hdr)

        # Create the information in the asn file
        
        asn_mem_names = np.array([''])
        asn_mem_types = np.array(['EXP-DTH', 'PROD-DTH'])
        asn_mem_prsnt = np.array([1,1,1,])

        hdu_data = fits.BinTableHDU().from_columns([fits.Column(name='MEMNAME', format='14A', array=asn_mem_names), 
                    fits.Column(name='MEMTYPE', format='14A', array=asn_mem_types), 
                    fits.Column(name='MEMPRSNT', format='L', array=asn_mem_prsnt)])

        # Create the final asn file
        hdu = fits.HDUList([primary_hdu, hdu_data])
        hdu.writeto(asn_filename)


    

        pass

    def run_reduction():
        pass

def get_flat(file_name):
    ''' Will check if user has proper reference file directories 
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
    
    os.environ['jref'] = '~/jref/'
    if not os.path.exists('jref'):
        os.mkdir('jref')

    base_url = 'https://hst-crds.stsci.edu/unchecked_get/references/hst/'
    
    with fits.open(file_name) as fitsfile:
        reffile_name = hdu[0].header['PFLTFILE'].replace('$', '/')
        if not os.path.exists(reffile_name):
            urlretrieve(base_url + os.path.basename(reffile_name), reffile_name)

    return reffile_name

def main():
    ''' Main function of reduce_dash. 

    Parameters
    ----------
     : 
        

    Outputs
    ----------
     : 
         

    '''


