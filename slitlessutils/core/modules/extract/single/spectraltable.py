import numpy as np
from astropy.stats import sigma_clip

from ....utilities import indices


class SpectralTable(dict):
    '''
    Class to hold measurements of a single-orient spectrum.

    inherits from `dict`

    Parameters
    ----------
    segid : int
       The segmentation ID for the source in question

    '''

    def __init__(self, segid):
        self.segid = segid
        self['wave'] = []  # wavelength vector
        self['dwav'] = []  # width of pixel in wavelength
        self['flam'] = []  # spectrum in CGS
        self['func'] = []  # uncertainty in spectrum in CGS
        self['flux'] = []  # count rate in e-/s
        self['cont'] = []  # contamination in CGS
        self['npix'] = []  # number of extracted pixels
        self.nrows = 0
        self.ncols = len(self.keys())

    def extend(self, *args):
        '''
        Method to grow the table

        Parameters
        ----------
        args : `tuple`
           Inputs of same length as ncols

        '''
        if self.ncols == len(args):
            n = len(args[0])

            if all(n == len(a) for a in args[1:]):

                for k, a in zip(self.keys(), args):
                    self[k].extend(a)
                self.nrows += n

    def append(self, *args):
        '''
        Method to append to the table

        Parameters
        ----------
        args : `tuple`
           Inputs of same length as ncols

        '''
        if self.ncols == len(args):
            for k, a in zip(self.keys(), args):
                self[k].append(a)
            self.nrows += 1

    def clear(self):
        '''
        Method to clear the table
        '''
        self.nrows = 0
        super().clear()

    def get(self, k):
        '''
        Method to get a column

        Parameters
        ----------
        k : `str`
           The name of the column
        '''

        if k in self:
            return np.asarray(self[k])
        else:
            raise AttributeError(f'Invalid keyword: {k}')

    def __len__(self):
        '''
        Return the number of rows
        '''
        return self.nrows

    def __add__(self, this):
        '''
        Method to concatenate two tables

        Parameters
        ----------
        this : `SpectralTable`
           The other table to include
        '''

        if this.segid == self.segid:
            new = SpectralTable(self.segid)

            for k in self.keys():
                new[k] = self[k] + this[k]

            return new

    def __iadd__(self, this):
        '''
        Method to concatenate two tables

        Parameters
        ----------
        this : `SpectralTable`
           The other table to include
        '''

        if this.segid == self.segid:
            for k in self.keys():
                self[k].extend(this[k])

            return self

    def __str__(self):
        return f'SpectralTable for object #{self.segid}'

    def __getattr__(self, k):
        '''
        Method to get a column

        Parameters
        ----------
        k : `str`
           The name of the column
        '''

        if k in self:
            return np.asarray(self[k])
        else:
            raise AttributeError(f'Invalid keyword: {k}')

    def combine(self, pars, clipper=None):
        '''
        Method to stack spectra

        Parameters
        ----------

        pars : `Disperser`
           A light-weight object that describes the extraction wavelengths

        clipper : `astropy.stats.SigmaClip`
           A method to sigma clip before combining.  Default is `None`

        Returns
        -------
        spectrum : `SpectralTable`
           A new `SpectralTable` object of the combined spectrum
        '''

        spectrum = SpectralTable(self.segid)

        # dump the results
        wave = self.get('wave')
        flam = self.get('flam')
        func = self.get('func')
        flux = self.get('flux')
        cont = self.get('cont')
        npix = self.get('npix')

        # bin to the output grid
        outw = pars.wavelengths()
        lamb = pars.indices(wave)
        ri = indices.reverse(lamb)

        for ind, g in ri.items():
            # extract the relevant data
            g = g[0]
            thisflam = flam[g]
            thisfunc = func[g]
            thisflux = flux[g]
            thiscont = cont[g]
            thisnpix = npix[g]

            # get the valid measurements
            g = np.isfinite(thisflam)
            thisflam = thisflam[g]
            thisfunc = thisfunc[g]
            thisflux = thisflux[g]
            thiscont = thiscont[g]
            thisnpix = thisnpix[g]

            # sigma clip and remove bad data
            if clipper:
                clipped = clipper(thisflam)
                g = np.logical_not(clipped.mask)
                thisflam = thisflam[g]
                thisfunc = thisfunc[g]
                thisflux = thisflux[g]
                thiscont = thiscont[g]
                thisnpix = thisnpix[g]

            # compute weights
            w = 1. / thisfunc**2

            # compute the weighted moments
            wsum = np.nansum(w)
            if wsum > 0:
                flamsum = np.nansum(w * thisflam) / wsum
                fluxsum = np.nansum(w * thisflux) / wsum
                flamvar = 1. / np.sqrt(wsum)
            else:
                flamsum = np.inf
                fluxsum = np.inf
                flamvar = np.inf

            numb = np.sum(thisnpix)
            contam = 0.0

            # save the results
            spectrum.append(outw[ind], pars.dwave, flamsum,
                            np.sqrt(flamvar), fluxsum, contam, numb)

        return spectrum

    def to_csv(self, filename=None):
        '''
        Method to write to a csv file

        Parameters
        ----------
        filename : `str`
            Name of the output file.  If set to `None`, then will create
            filename of the form `f'{segid}.csv'`.  Default: None

        '''

        # separator for the output file
        sep = ' , '

        # default output name
        if filename is None:
            filename = f'{self.segid}.csv'

        # open output file
        with open(filename, 'w') as fp:
            print('# ' + sep.join(self.keys()), file=fp)
            for args in zip(*self.values()):
                print(*args, file=fp, sep=sep)


if __name__ == '__main__':
    x = SpectralTable(1)
    x.extend(np.arange(10), np.arange(10),
             np.arange(10), np.arange(10),
             np.arange(10), np.arange(10))
    print(x.nrows)
    print(x.combine(1).ncols)
