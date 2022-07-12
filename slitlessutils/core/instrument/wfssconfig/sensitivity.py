from astropy.io import fits
import numpy as np


from ...photometry import Band

class Sensitivity(Band):
    def __init__(self,sensfile,senslimit=1e10):
        self.sensfile=sensfile
        self.senslimit=senslimit

        # read the file
        data,header=fits.getdata(self.sensfile,1,header=True)


        g=np.where(data['SENSITIVITY'] > self.senslimit)

        Band.__init__(self,data['WAVELENGTH'],data['SENSITIVITY'],
            where=g,unit=header.get('TUNIT1',''))


    def __str__(self):
        return f'Sensitivity curve: {self.sensfile}'
