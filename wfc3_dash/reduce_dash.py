from astropy.io import fits

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
        
    def run_reduction():
        
        pass