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

from astropy.io import fits

FLAT_F160W = fits.open(''))[1].data
BP_MASK = fits.open('')[0].data
PATH = 'data/'


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
        
    
    def split_ima(self):

        pass
        
    def subtract_background_reads():
        
        pass
        
    def fix_cosmic_rays():
        
        pass
        
    def align_reads():
        
        pass
        
    def coadd_reads():
        
        pass
        
    def run_reduction():
        
        pass