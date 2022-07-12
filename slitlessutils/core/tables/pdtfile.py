import os
import numpy as np
from shapely.geometry import Polygon
from shapely.ops import unary_union




from .hdf5file import HDF5File
from . import attributes
from .opt import OPT
from .pdt import PDT
from .rdt import RDT

from ..utilities import indices



class PDTFile(HDF5File):
    TTYPE='pdt'

    def __init__(self,wfssfile,**kwargs):
        HDF5File.__init__(self,wfssfile,**kwargs)


        
    def load_polygon(self,source):
        ''' load a `Source` as a Shapely polygon '''

        # create empty lists and fill them from a PDT
        x,y,l=[],[],[]
        for x,y in source.pixels():
            pdt=self.load_pdt(x,y)
            x.extend(pdt['x'])
            y.extend(pdt['y'])
            l.extend(pdt['l'])

        # make them arrays to make for better integration with shapely
        x=np.array(x)
        y=np.array(y)
        l=np.array(l)


        # collect based on wavelength
        ri=indices.reverse(l)
        polys=[Polygon(list(zip(x[g],y[g]))) for ll,g in ri.items()]

        # group the polygons
        poly = unary_union(polys)

        return poly



    def load_omt(self,source):
        ''' load all the pixels in a `Source` as a mask '''

        # build an empty table and fill it pixel by pixel
        omt=OMT(source)
        for x,y in source.pixels():
            # load the order
            hd=self.h5order[omt.name]

            # copy the data into the columns of the output PDT
            data=hd[()]
            omt.extend(data['x'],data['y'])

        # only keep unique WFS image pixels
        omt.uniqify()

        # load and copy over any attributes (think header keywords)
        for attr in hd.attrs:
            if attr not in ('x','y'):
                omt.attrs[attr]=attributes.load(hd,attr)

        return omt

    def load_odt(self,source):

        odt=ODT(source)
        for x,y in source.pixels:
            xd,yd=source.image_coordinate(x,y,dtype=int)
            pdt=self.load(xd,yd)
            odt.append(pdt)

        odt.decimate()


        for k,v in pdt.attrs.items():
            if k not in ('x','y'):
                odt.attrs[k]=v

        return odt


    def load_opt(self,source):
        ''' load an Object Profile Table (OPT) '''
        
        opt=OPT(source.segid)
        for x,y in source.pixels():
            xd,yd=source.image_coordinates(x,y,dtype=int)
            pdt=self.load_pdt(xd,yd)
            opt.append(pdt)

        opt.decimate()

        for k,v in pdt.attrs.items():
            if k not in ('x','y'):
                opt.attrs[k]=v
                
        
        return opt
        


    
    def load_rdts(self,source):
        ''' load the Region Dispersion Tables (RDTs) '''
        rdts={}
        for regid,region in enumerate(source):
            rdt=RDT(source,regid)
            
            for x,y in region.pixels():
                pdt=self.load_pdt(x,y)
                if len(pdt)==0:
                    LOGGER.debug("MAJOR ERROR IN READING REGION")
                rdt.append(pdt)
            rdt.decimate()

            for k,v in pdt.attrs.items():
                if k not in ('x','y'):
                    rdt.attrs[k]=v

            rdts[regid]=rdt
        return rdts

    def load_pdts(self,source,flatten=False):
        ''' load multiple Pixel Dispersion Tables (PDTs) at once '''

        pdts={}
        if flatten:
            for x,y in source.pixels():
                xyd=source.image_coordinates(x,y,dtype=int)
                pdts[xyd]=self.load_pdt(*xyd)
        else:
            for regid,region in source.items():
                pdts[regid]={}
                for x,y in region.pixels():
                    xyd=source.image_coordinates(x,y,dtype=int)
                    pdts[regid][xyd]=self.load_pdt(*xyd)        
                            
        return pdts

    def load_pdt(self,x,y):
        ''' load a Pixel Dispersion Table (PDT) '''

        pdt = PDT(x,y)

        # check if this order is present in the PDTFile
        if pdt.name in self.h5order:
            # load the order
            hd=self.h5order[pdt.name]

            # load and copy over any attributes (think header keywords)
            for attr in hd.attrs:
                pdt.attrs[attr]=attributes.load(hd,attr)

            # copy the data into the columns of the output PDT
            data=hd[()]
            for column in pdt.COLUMNS:
                pdt[column]=data[column]

        return pdt
