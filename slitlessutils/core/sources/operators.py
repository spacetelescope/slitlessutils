import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage as ndi
from skimage import morphology as morph

from slitlessutils.logger import LOGGER

from ..utilities import headers


class Operator:
    """
    Base class to implement image operators
    """

    def __str__(self):
        return self.__str

    def __setattr__(self, k, v):
        if k in self.__pars:
            self.__pars[k] = (v, self.__pars[k][1])

    def __getattr__(self, k):
        if k == "_Operator__pars":
            return self.__pars
        if k in self.__pars:
            return self.__pars[k][0]

    def __and__(self, obj):
        new = OperatorCollection()
        if isinstance(self, list):
            new.extend(self)
        else:
            new.append(self)
        if isinstance(obj, list):
            new.extend(obj)
        else:
            new.append(obj)
        return new

    def update_header(self, hdr):
        """
        Method to update a fits header

        Parameters
        ----------
        hdr : `astropy.io.fits.Header`
            The fits header to update

        """

        before = None
        for k, v in self.__pars.items():
            hdr[k] = v
            if not before:
                before = k
        headers.add_stanza(hdr, str(self), before=before)


class OperatorCollection(Operator, list):
    """
    Class to chain together multiple operators
    """

    def __str__(self):
        s = 'Collection of morphological operators: '
        if len(self) == 0:
            return s + '<EMPTY>'
        else:
            return s + '\n'.join(str(s) for s in self)

    def __call__(self, *args):
        """
        Method to call the operator

        Parameters
        ----------
        args : tuple
            All operators take as input the direct image, segmentation
            image, and direct image header, in that order.

        Returns
        -------
        args : tuple
            All operators returned the modified direct image, modified
            segmentation image, and modified direct image header, in
            that order.  Therefore, the outputs of one operator can be
            inputs to a new operator.

        Notes
        -----
        The inputs/outputs could be set in this method directly, but I
        chose to code it this way so that the operators can change
        underneath this driver code.
        """

        for op in self:
            args = op(*args)
        return args

    def update_header(self, hdr):
        """
        Method to update a fits header

        Parameters
        ----------
        hdr : `astropy.io.fits.Header`
            The fits header to update
        """
        for op in self:
            op.update_header(hdr)


class Erode(Operator):
    """
    Class to implement the erosion operator using `skimage.morphology`.

    This will not affect the direct image, as it only erodes the segmentation
    map.

    Notes
    -----
    There are three erosion kernels: 'disk','square',' and 'diamond',
    that can be affected by the class variable 'erodefun'.  The spatial
    scale of the dilation can be changed with 'eroderad'

    Examples
    --------
    >>> e=Erode()
    >>> e.eroderad=4  # doctest: +SKIP
    >>> e.erodefun='diamond'  # doctest: +SKIP
    >>> img,seg,hdr = e(img,seg,hdr)  # doctest: +SKIP
    """

    __str = 'Erosion Operator'
    __pars = {'eroderad': (2, 'size of erosion kernel'),
              'erodefun': ('disk', 'shape of erosion kernel')}
    __shapes = {'disk': morph.disk, 'square': morph.square, 'diamond': morph.diamond}

    def __call__(self, img, seg, hdr):
        """
        Method to call the erosion

        Parameters
        ----------
        img : `np.ndarray`
            The direct image.

        seg : `np.ndarray`
            The segmentation image.

        hdr : `astropy.io.fits.Header()`
            The direct image header

        Returns
        -------
        newimg : `np.ndarray`
            The direct image.

        newseg : `np.ndarray`
            The eroded segmentation image.

        newhdr : `astropy.io.fits.Header()`
            The updated direct image header
        """

        if self.eroderad > 0 and self.erodefun in self.__shapes:
            fun = self.__shapes[self.erodefun.lower()]
            footprint = fun(self.eroderad)

            # find the pixels with this segid
            gpx = seg == hdr['SEGID']

            # erode the image
            newseg = morph.binary_erosion(gpx, footprint=footprint).astype(np.int)

            # find pixels
            b = np.where((seg != 0) & (seg != hdr['SEGID']))
            g = np.where(newseg)
            newseg[g] = seg[g]
            newseg[b] = seg[b]

            # axes[0].imshow(seg,origin='lower')
            # axes[1].imshow(newseg,origin='lower')
            # axes[2].imshow(seg-newseg,origin='lower')
            # plt.show()

            # record what was done
            newhdr = hdr.copy()
            self.update_header(newhdr)

            return img, newseg, newhdr
        else:
            LOGGER.warning('Unable to erode.')
            return img, seg, hdr


class Dilate(Operator):
    """
    Class to implement the dilation operator using `skimage.morphology`.

    This will not affect the direct image, as it only dilates the
    segmentation map.

    Notes
    -----
    There are three dilation kernels: 'disk','square',' and 'diamond',
    that can be affected by the class variable 'dilatfun'.  The spatial
    scale of the dilation can be changed with 'dilatrad'

    Examples
    --------
    >>> d=Dilate()
    >>> d.dilatrad=4  # doctest: +SKIP
    >>> d.dilatfun='diamond'  # doctest: +SKIP
    >>> img,seg,hdr = d(img,seg,hdr)  # doctest: +SKIP

    """

    __str = 'Dilation Operator'
    __pars = {'dilatrad': (3, 'size of dilation kernel'),
              'dilatfun': ('square', 'shape of dilation kernel')}
    __shapes = {'disk': morph.disk, 'square': morph.square, 'diamond': morph.diamond}

    def __call__(self, img, seg, hdr):
        """
        Method to call the erosion

        Parameters
        ----------
        img : `np.ndarray`
            The direct image.

        seg : `np.ndarray`
            The segmentation image.

        hdr : `astropy.io.fits.Header()`
            The direct image header

        Returns
        -------
        newimg : `np.ndarray`
            The direct image.

        newseg : `np.ndarray`
            The dilated segmentation image.

        newhdr : `astropy.io.fits.Header()`
            The updated direct image header
        """

        if self.dilatrad > 0 and self.dilatfun in self.__shapes:
            fun = self.__shapes[self.dilatfun.lower()]
            footprint = fun(self.dilatrad)

            # find the pixels with this segid
            gpx = seg == hdr['SEGID']

            # erode the image
            newseg = morph.binary_dilation(gpx, footprint=footprint).astype(np.int)

            # find pixels
            b = np.where((seg != 0) & (seg != hdr['SEGID']))
            g = np.where(newseg)
            newseg[g] = hdr['SEGID']
            newseg[b] = seg[b]

            # record what was done
            newhdr = hdr.copy()
            self.update_header(newhdr)

            return img, newseg, newhdr
        else:
            LOGGER.warning('Unable to dilate.')
            return img, seg, hdr


class Smooth(Operator):
    """
    Class to implement the image smoothing

    This will not affect the segmentation image, as it only smoothes
    the direct image.

    Notes
    -----
    There are two options for smoothing: median and gaussian.  They can
    changed by changing the class variable 'smthfunc'.  The spatial scale
    can be changed with 'smthsize'.

    Examples
    --------
    >>> s=Smooth()
    >>> s.smthsize=4  # doctest: +SKIP
    >>> s.smthfunc='median'  # doctest: +SKIP
    >>> img,seg,hdr = s(img,seg,hdr)  # doctest: +SKIP

    """

    __str = 'Smoothing Operator'
    __pars = {'smthsize': (2, 'Scale of smoothing operation'),
              'smthfunc': ('guassian', 'Smoothing function')}

    __funcs = {'gaussian': ndi.gaussian_filter,
               'median': ndi.median_filter}

    def __call__(self, img, seg, hdr):
        """
        Method to call the smoothing

        Parameters
        ----------
        img : `np.ndarray`
            The direct image.

        seg : `np.ndarray`
            The segmentation image.

        hdr : `astropy.io.fits.Header()`
            The direct image header

        Returns
        -------
        newimg : `np.ndarray`
            The smoothed direct image.

        newseg : `np.ndarray`
            The segmentation image.

        newhdr : `astropy.io.fits.Header()`
            The updated direct image header
        """

        if self.smthsize > 0 and self.smthfunc in self.__funcs:
            filt = self.__funcs[self.smthfunc.lower()]

            bad = seg != hdr['SEGID']

            V = img.copy()
            V[bad] = 0
            VV = filt(V, self.smthsize)

            W = np.ones_like(img)
            W[bad] = 0
            WW = filt(W, self.smthsize)
            new = VV / np.maximum(WW, 1e-5)

            # fig,axes=plt.subplots(2,2)
            # axes[0][0].imshow(img,origin='lower')
            # axes[0][1].imshow(new,origin='lower')
            # axes[1][0].imshow(VV,origin='lower')
            # axes[1][1].imshow(WW,origin='lower')
            # plt.show()

            newhdr = hdr.copy()
            self.update_header(newhdr)

            return new, seg, newhdr
        else:
            LOGGER.warning('Unable to smooth.')
            return img, seg, hdr


class Rebin(Operator):
    """
    class to rebin the direct image and segmentation map

    Notes
    -----
    There are three rebinning functions: 'mean', 'median', and 'sum'
    that can be affected by the class variable 'binfun'.  The spatial
    factor of the rebinning is set by 'binfac'.

    Examples
    --------
    >>> r=Rebin()
    >>> r.binfac=2  # doctest: +SKIP
    >>> r.binfun='median'  # doctest: +SKIP
    >>> img,seg,hdr = r(img,seg,hdr)  # doctest: +SKIP

    """

    __str = 'Rebinning Operator'
    __pars = {'binfac': (2, 'Rebinning factor'),
              'binfun': ('mean', 'Rebinning function')}
    __funcs = {'mean': ndi.mean, 'median': ndi.median, 'sum': ndi.sum_labels}

    CDMAT = ("CD1_1", "CD1_2", 'CD2_1', 'CD2_2')
    # PCMAT=("PC1_1","PC1_2",'PC2_1','PC2_2')
    # PCVEC=("PC1_1","PC2_2")
    # CDELT=('CDELT1','CDELT2')

    def __call__(self, img, seg, hdr):
        """
        Method to rebin the image and segmentation map

        Parameters
        ----------
        img : `np.ndarray`
            The direct image.

        seg : `np.ndarray`
            The segmentation image.

        hdr : `astropy.io.fits.Header()`
            The direct image header

        Returns
        -------
        newimg : `np.ndarray`
            The rebinned direct image.

        newseg : `np.ndarray`
            The rebinned segmentation image.

        newhdr : `astropy.io.fits.Header()`
            The updated direct image header
        """

        binfact = self.binfact
        if binfact > 1 and self.binfun in self.__funcs:

            dim = img.shape
            new = [(dim[0] + 1) // binfact, (dim[1] + 1) // binfact]
            bins = np.arange(dim[0] * dim[1]).reshape(dim) + 1
            ones = np.ones((binfact, binfact), dtype=np.int)
            bins = np.kron(bins, ones)[:dim[0], :dim[1]]
            b = np.unique(bins)

            func = self.__funcs[self.binfun.lower()]
            newimg = func(img, labels=bins, index=b).reshape(new)

            newhdr = hdr.copy()
            newhdr['NAXIS1'] = newimg.shape[1]
            newhdr['NAXIS2'] = newimg.shape[0]

            newhdr['CRPIX1'] = (newhdr['CRPIX1'] + 0.5) / binfact + 0.5
            newhdr['CRPIX2'] = (newhdr['CRPIX2'] + 0.5) / binfact + 0.5

            if all(k in newhdr for k in self.CDMAT):
                newhdr['CD1_1'] *= binfact
                newhdr['CD2_1'] *= binfact
                newhdr['CD1_2'] *= binfact
                newhdr['CD2_2'] *= binfact
            else:  # if all(k in newhdr for k in self.PCMAT):
                newhdr['CDELT1'] *= binfact
                newhdr['CDELT2'] *= binfact

            LOGGER.debug('REBINNING SEGMAP IS PROBABLY NOT RIGHT')
            g = np.where(seg != newhdr['SEGID'])    # pixels not the source
            newseg = seg.copy()
            newseg[g] = 0
            newseg = ndi.maximum(newseg, labels=bins, index=b).reshape(new)
            tmp = ndi.maximum(seg, labels=bins, index=b).reshape(new)
            g = np.where(newseg == 0)[0]
            newseg[g] = tmp[g]

            fig, axes = plt.subplots(2, 2)
            axes[0][0].imshow(img, origin='lower')
            axes[0][1].imshow(seg, origin='lower')
            axes[1][0].imshow(newseg, origin='lower')
            # axes[1][1].imshow(WW,origin='lower')
            plt.show()

            self.update_header(newhdr)

            return newimg, newseg, newhdr
        else:
            LOGGER.warning('Unable to rebin.')
            return img, seg, hdr


if __name__ == '__main__':
    e = Erode()
    s = Smooth()
    r = Rebin()
    a = e & s & r

    q = e(1, 2, 3)
