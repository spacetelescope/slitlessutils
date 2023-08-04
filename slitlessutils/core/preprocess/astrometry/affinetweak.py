from astropy.io import fits
import numpy as np
import os


class AffineTweak(dict):
    def __init__(self, reffile, key0='A', key1=''):

        self.reffile = reffile
        self.key0 = key0
        self.key1 = key1

        with fits.open(self.reffile, mode='readonly') as hdul:
            for hdu in hdul:
                extname = hdu.header.get('EXTNAME')
                extver = hdu.header.get('EXTVER')

                if extname == 'SCI':
                    cd0 = self.get_cd(hdu.header, self.key0)
                    crval0 = self.get_crval(hdu.header, self.key0)

                    cd1 = self.get_cd(hdu.header, self.key1)
                    crval1 = self.get_crval(hdu.header, self.key1)

                    ext = (extname, extver)
                    self[ext] = {}
                    self[ext]['A'] = np.dot(cd1, np.linalg.inv(cd0))
                    self[ext]['d'] = crval1-crval0
                    self[ext]['wcsname0'] = hdu.header[f'WCSNAME{self.key0}']
                    self[ext]['wcsname1'] = hdu.header[f'WCSNAME{self.key1}']

    def __call__(self, datfile, inplace=False, newfile=None):

        mode = 'update' if inplace else 'readonly'
        with fits.open(datfile, mode=mode) as hdul:
            for ext, data in self.items():

                cd = self.get_cd(hdul[ext].header, self.key0)
                crval = self.get_crval(hdul[ext].header, self.key0)

                cd = np.dot(data['A'], cd)
                crval = crval+data['d']

                self.set_cd(hdul[ext].header, cd, self.key1)
                self.set_crval(hdul[ext].header, crval, self.key1)

            hdul[0].header.add_history('Affine tweaked astrometry')

            if inplace:
                outfile = datfile
            else:
                if newfile is None:
                    base = os.path.splitext(os.path.basename(datfile))[0]
                    newfile = f'{base}_twcs.fits'
                hdul.writeto(newfile, overwrite=True)
                outfile = newfile
        return outfile

    @staticmethod
    def get_cd(hdr, key):
        cd = np.array([[hdr[f'CD1_1{key}'], hdr[f'CD1_2{key}']],
                       [hdr[f'CD2_1{key}'], hdr[f'CD2_2{key}']]],
                      dtype=float)
        return cd

    @staticmethod
    def set_cd(hdr, cd, key):
        hdr[f'CD1_1{key}'] = cd[0, 0]
        hdr[f'CD1_2{key}'] = cd[0, 1]
        hdr[f'CD2_1{key}'] = cd[1, 0]
        hdr[f'CD2_2{key}'] = cd[1, 1]

    @staticmethod
    def get_crval(hdr, key):
        crval = np.array([hdr[f'CRVAL1{key}'], hdr[f'CRVAL2{key}']],
                         dtype=float)
        return crval

    @staticmethod
    def set_crval(hdr, crval, key):
        hdr[f'CRVAL1{key}'] = crval[0]
        hdr[f'CRVAL2{key}'] = crval[1]

    # def __str__(self):
    #    print(
