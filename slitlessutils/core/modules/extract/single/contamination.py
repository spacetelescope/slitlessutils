from astropy.io import fits
import numpy as np

from ....utilities import indices, headers
from .....logger import LOGGER


class ContOrders:
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

        if isinstance(self.cntorders, str) and self.cntorders.lower() == 'all':
            yield from h5.h5detector.keys()

        elif isinstance(self.cntorders, (tuple, list)):
            yield from self.cntorders
        else:
            pass

    def update_header(self, hdr):
        hdr.set('CONORDS', value=str(self.cntorders), comment='orders masked')


class Contamination(dict):

    def __init__(self, cntorders):
        self.cntorders = ContOrders(cntorders)
        self.threshold = 0.0

    def update_header(self, hdr):
        hdr.set('CONTAM', value=True, comment='was contamination applied')
        self.cntorders.update_header(hdr)
        hdr.set('CONTHRSH', value=self.threshold, comment='threshold contamination (e-/s)')
        headers.add_stanza(hdr, 'Contamination Settings', before='CONTAM')

    def make_model(self, sources, h5):
        LOGGER.info(f'Building contamination model for {h5.dataset}')

        # load all the polygons
        for ordname in self.cntorders(h5):
            h5.load_order(ordname)
            for segid, source in sources.items():
                poly = h5.load_polygon(source)
                if poly:
                    self[(segid, ordname)] = poly

    # def plotpolygons(self,segid,ordname):
    #    for (seg,ord),poly in self.items():
    #        if ordname==ord and seg==segid:
    #
    #
    # def plot(self):
    #    for (segid,ordname),v in self.items():
    #        if ordname=='+1':
    #            plt.plot(*v.exterior.xy,label=str(segid))
    #    plt.legend()
    #    plt.show()

    def __call__(self, segid, ordname, h5, sources, detdata, flatfield):

        # grab the polygon
        key = (segid, ordname)
        poly = self[key]

        # get teh bounds
        bounds = poly.bounds
        x0 = np.maximum(int(np.floor(bounds[0]))-1, 0)
        y0 = np.maximum(int(np.floor(bounds[1]))-1, 0)
        x1 = np.minimum(int(np.ceil(bounds[2]))+1, detdata.naxis[0]-1)
        y1 = np.minimum(int(np.ceil(bounds[3]))+1, detdata.naxis[1]-1)

        # create an image of all zeros
        dx = x1-x0+1
        dy = y1-y0+1
        img = np.zeros((dy, dx), dtype=float)

        # get a header and fill it
        hdr = detdata.wcs.to_header(relax=True)

        # update the header with keywords for the
        hdr.insert(0, ("XTENSION", 'IMAGE', 'Image extension'))
        hdr.insert(1, ('BITPIX', -64, 'array data type'))
        hdr.insert(2, ('NAXIS', 2, 'number of array dimensions'))
        hdr.insert(3, ('NAXIS1', dx))
        hdr.insert(4, ('NAXIS2', dy))
        hdr.insert(5, ('PCOUNT', 0, 'number of parameters'))
        hdr.insert(6, ('GCOUNT', 1, 'number of groups'))

        # update CRPIX for the cutout
        hdr['CRPIX1'] -= x0
        hdr['CRPIX2'] -= y0

        # give the LTV keywords for physical and image coordiantes
        hdr['LTV1'] = (-x0, 'x-offset between phys. and image')
        hdr['LTV2'] = (-y0, 'y-offset between phys. and image')

        # add an extname and a few things
        hdr['EXTNAME'] = f'{detdata.name},{segid},{ordname}'
        hdr['DETNAME'] = (detdata.name, 'Detector name')
        hdr['SEGID'] = (segid, 'Segmentation ID for this cont model')
        hdr['ORDNAME'] = (ordname, 'Order for cont model')

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

                            # sum over unique pixels for the weights and
                            # weighted-average for the wavelength
                            vv, yy, xx = indices.decimate(val, yg, xg,
                                                          dims=detdata.shape)
                            ww, _y, _x = indices.decimate(val*wav, yg, xg,
                                                          dims=detdata.shape)
                            ww /= vv

                            # get some correction factors
                            flat = flatfield(xx, yy, ww)
                            sens = detdata.config[ordname].sensitivity(ww)
                            flam = region.sed(ww, fnu=False)
                            vv *= (flat*sens*flam)

                            # flat=flatfield(xg,yg,wav)
                            # sens=detdata.config[ordname].sensitivity(wav)
                            # flam=region.sed(wav,fnu=False)
                            # val*=(sens*flat*flam)
                            #
                            # vv,yy,xx=indices.decimate(val,yg,xg,
                            #                          dims=detdata.shape)

                            # only work with pixels that are inside the cutout
                            g = np.where((x0 <= xx) & (xx < x1+1) &
                                         (y0 <= yy) & (yy < y1+1))[0]

                            if g.size > 0:
                                xx = xx[g]
                                yy = yy[g]
                                vv = vv[g]

                                # update for the relative pixel area
                                vv *= detdata.relative_pixelarea(xx, yy)

                                # apply a threshold to the levels
                                vv = np.where(vv > self.threshold, vv, 0.0)

                                # sum in to the image
                                img[yy-y0, xx-x0] += vv

        # return the data as fits.ImageHDU so we can do stuff later
        return fits.ImageHDU(data=img, header=hdr)
