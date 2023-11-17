import numpy as np
from astropy.io import fits

from .....logger import LOGGER
from ....utilities import headers, indices


class ContOrders:
    """
    Class to contain info for the spectral orders considered contaminants

    Parameters
    ----------
    cntorders : list/tuple
        the orders to consider contaminants
    """

    def __init__(self, cntorders):
        self.cntorders = cntorders

    def __str__(self):
        if isinstance(self.cntorders, str) and self.cntorders.lower() == 'all':
            conorders = '<all>'
        elif isinstance(self.cntorders, (tuple, list)):
            if len(self.cntorders) == 0:
                conorders = '<all>'
            else:
                conorders = ','.join(self.cntorders)
        else:
            conorders = '<None>'

        return conorders

    def __call__(self, h5):
        """
        A generator that iterates over the orders and returns HDF5 dataset

        Parameters
        ----------
        h5 : an HDF5 object
            The base of the HDF5 filesystem

        Returns
        -------
        The HDF5 dataset
        """

        if isinstance(self.cntorders, str) and self.cntorders.lower() == 'all':
            yield from h5.h5detector.keys()

        elif isinstance(self.cntorders, (tuple, list)):
            yield from self.cntorders
        else:
            pass

    def update_header(self, hdr):
        hdr.set('CONORDS', value=str(self.cntorders), comment='orders masked')


class Contamination(dict):
    """
    Class to compute 1d contamination

    Parameters
    ----------
    cntorders : list/tuple or str
        The orders to consider as contamination.  If "all" then all
        orders are masked.

    """

    def __init__(self, cntorders):
        self.cntorders = ContOrders(cntorders)
        self.threshold = 0.0

    def update_header(self, hdr):
        """
        Method to update a fits header

        Parameters
        ----------
        hdr : `astropy.io.fits.Header`
            The header to update
        """

        hdr.set('CONTAM', value=True, comment='was contamination applied')
        self.cntorders.update_header(hdr)
        hdr.set('CONTHRSH', value=self.threshold, comment='threshold contamination (e-/s)')
        headers.add_stanza(hdr, 'Contamination Settings', before='CONTAM')

    def make_model(self, sources, h5):
        """
        Method to initialize a contamination model

        Parameters
        ----------
        sources : `su.sources.SourceCollection`
            The collection of sources to model

        h5 : HDF5 object
            The HDF5 dataset that contains the data for this detector

        Returns
        -------
        None
        """

        LOGGER.info(f'Building contamination model for {h5.dataset}{h5.h5detector.name}')

        # load all the polygons
        for ordname in self.cntorders(h5):
            h5.load_order(ordname)
            for segid, source in sources.items():
                poly = h5.load_polygon(source)
                if poly:
                    self[(segid, ordname)] = poly

    def __call__(self, segid, ordname, h5, sources, detdata, flatfield, bbx=None, bby=None):
        """
        Method to compute the contamination model

        Parameters
        ----------
        segid : int
            The segmentation ID to compute contamination model

        ordname : str
            The spectral order for which the model is computed

        h5 : HDF5 object
            The HDF5 dataset were the forward model is stored

        sources : `su.source.SourceCollection`
            The collection of spectral sources

        detdata : `su.wfss.data.DetectorData`
            The detector data for this contamination model

        flatfield : `su.wfss.config.FlatField`
            The flat field to use.

        bbx : 2-tuple or None
            The bounding box for the x.  If None, then will use full frame.  Default is None.

        bby : 2-tuple or None
            The bounding box for the y.  If None, then will use full frame.  Default is None.

        Returns
        -------
        dat : `astropy.io.fits.ImageHDU`
            The image and header for the contamination

        """

        # grab the polygon
        key = (segid, ordname)
        poly = self[key]

        if bbx is None:
            bbx = (0, detdata.naxis[0]-1)
        nx = bbx[1]-bbx[0]+1

        if bby is None:
            bby = (0, detdata.naxis[1]-1)
        ny = bby[1]-bby[0]+1

        # an image of the contamination
        img = np.zeros((ny, nx), dtype=float)

        # process each polygon
        for k, ply in self.items():

            # ignore if it's the polygon in question
            if k != key:

                # find out if it overlaps
                dp = poly.intersection(ply)
                if not dp.is_empty:

                    # load the data for this order
                    h5.load_order(ordname)

                    # process each spectral region
                    for regid, region in sources[k[0]].items():
                        for x, y in region.pixels(applyltv=True, dtype=int):
                            # load the data
                            pdt = h5.load_pdt(x, y)

                            xg = pdt.get('x')
                            yg = pdt.get('y')
                            val = pdt.get('val')
                            wav = pdt.wavelengths()

                            g = np.where(
                                (bbx[0] <= xg) & (
                                    xg < bbx[1]) & (
                                    bby[0] <= yg) & (
                                    yg < bby[1]))[0]
                            if g.size > 0:
                                xg = xg[g]
                                yg = yg[g]
                                val = val[g]
                                wav = wav[g]

                                # apply the calibrations for this pixel
                                flat = flatfield(xg, yg, wav)
                                sens = detdata.config[ordname].sensitivity(wav)
                                area = detdata.relative_pixelarea(xg, yg)
                                flam = region.sed(wav, fnu=False)
                                val *= (sens*area*flat*flam)

                                # sum over pixels
                                val, xx, yy = indices.decimate(val, xg, yg)
                                img[yy-bby[0], xx-bbx[0]] += val

        # return the data as fits.ImageHDU so we can do stuff later
        hdr = detdata.wcs.to_header(relax=True)

        # update the header with keywords for the
        hdr.insert(0, ("XTENSION", 'IMAGE', 'Image extension'))
        hdr.insert(1, ('BITPIX', -64, 'array data type'))
        hdr.insert(2, ('NAXIS', 2, 'number of array dimensions'))
        hdr.insert(3, ('NAXIS1', img.shape[1]))
        hdr.insert(4, ('NAXIS2', img.shape[0]))
        hdr.insert(5, ('PCOUNT', 0, 'number of parameters'))
        hdr.insert(6, ('GCOUNT', 1, 'number of groups'))

        # update CRPIX for the cutout
        hdr['CRPIX1'] -= bbx[0]
        hdr['CRPIX2'] -= bby[0]

        # give the LTV keywords for physical and image coordinates
        hdr['LTV1'] = (-bbx[0], 'x-offset between phys. and image')
        hdr['LTV2'] = (-bby[0], 'y-offset between phys. and image')

        # add an extname and a few things
        hdr['EXTNAME'] = f'{detdata.name},{segid},{ordname}'
        hdr['DETNAME'] = (detdata.name, 'Detector name')
        hdr['SEGID'] = (segid, 'Segmentation ID for this cont model')
        hdr['ORDNAME'] = (ordname, 'Order for cont model')

        return fits.ImageHDU(data=img, header=hdr)
