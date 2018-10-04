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
import numpy as np
from urllib.request import urlretrieve 

BP_MASK = fits.open('')[0].data
PATH = 'data_download/mastDownload/HST/'

class DashData:
    
    def __init__(self,file_name):
        
        if '_ima.fits' not in file_name:  #probably the header would be a better check?
            raise Exception('Input needs to be an IMA file.')
        else:
            self.file_name = file_name
            try:
                self.ima_file = fits.open(self.file_name)
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

            var = 2*self.readnoise_2D + science_data*FLAT    
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

        
    def subtract_background_reads():
        
        pass
        
    def fix_cosmic_rays():
        
        pass
        
    def align_reads():
        
        pass
        
    def coadd_reads():
        
        pass

    def make_pointing_asn(self):
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


