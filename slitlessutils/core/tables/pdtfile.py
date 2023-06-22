from shapely.geometry import Polygon

from . import attributes
from .hdf5file import HDF5File
from .opt import OPT
from .odt import ODT
from .pdt import PDT
from .rdt import RDT
from slitlessutils.logger import LOGGER


class PDTFile(HDF5File):
    """
    A class to hold a file of pixel-dispersion tables.

    inherits from `HDF5File`

    """

    # the type of the table. Ok, not really needed, but preparing for
    # some new files later in life?
    TTYPE = 'pdt'

    def __init__(self, wfssfile, **kwargs):
        """
        Initializer

        Parameters
        ----------
        wfssfile : `su.core.wfss.data.WFSS`
            The WFSS file to write/open/read a `PDTFile` for

        kwargs : dict, optional
            The optional keywords to pass to the `HDF5File` parent class
        """

        HDF5File.__init__(self, wfssfile, **kwargs)

    def load_polygon(self, source, **kwargs):
        """
        Method to load a shapely polygon from the PDT

        Parameters
        ----------
        source : `su.core.sources.Source`
            The source to load a polygon as

        kwargs : dict, optional
            Optional keywords passed to `compute_vertices()`

        Returns
        -------
        poly : `shapely.geometry.Polygon`
            The shapely polygon

        """

        # create empty lists and fill them from a PDT
        odt = self.load_odt(source)

        # group the polygons
        # poly = unary_union(polys)
        px, py = odt.compute_vertices(**kwargs)
        if len(px) > 0:
            poly = Polygon(list(zip(px, py)))
        else:
            poly = None

        return poly

    # def load_omt(self,source):
    #    ''' load all the pixels in a `Source` as a mask '''
    #
    #    # build an empty table and fill it pixel by pixel
    #    omt=OMT(source)
    #    for x,y in source.pixels():
    #        # load the order
    #        hd=self.h5order[omt.name]
    #
    #        # copy the data into the columns of the output PDT
    #        data=hd[()]
    #        omt.extend(data['x'],data['y'])
    #
    #    # only keep unique WFS image pixels
    #    omt.uniqify()
    #
    #    # load and copy over any attributes (think header keywords)
    #    for attr in hd.attrs:
    #        if attr not in ('x','y'):
    #            omt.attrs[attr]=attributes.load(hd,attr)
    #
    #    return omt

    def load_odt(self, source):
        """
        Method to load an object-dispersion table (ODT) from this `PDTFile`
        given a source

        Parameters
        ----------
        source : `su.core.sources.Source`
            The source to load

        Returns
        -------
        odt : `su.core.tables.ODT`
            The object-dispersion table
        """

        odt = ODT(source)
        for x, y in source.pixels():
            xd, yd = source.image_coordinates(x, y, dtype=int)
            pdt = self.load_pdt(xd, yd)
            odt.append(pdt)

        odt.decimate()

        for k, v in pdt.attrs.items():
            if k not in ('x', 'y'):
                odt.attrs[k] = v

        return odt

    def load_opt(self, source):
        """
        Method to load an object-profile table (OPT) from this `PDTFile`
        given a source

        Parameters
        ----------
        source : `su.core.sources.Source`
            The source to load

        Returns
        -------
        odt : `su.core.tables.OPT`
            The object-profile table
        """
        opt = OPT(source)
        for x, y in source.pixels():
            xd, yd = source.image_coordinates(x, y, dtype=int)
            pdt = self.load_pdt(xd, yd)
            opt.append(pdt)

        opt.decimate()

        for k, v in pdt.attrs.items():
            if k not in ('x', 'y'):
                opt.attrs[k] = v

        return opt

    def load_rdts(self, source):
        """
        Method to load a bunch of region-dispersion tables (RDTs) from
        this `PDTFile` given a source

        Parameters
        ----------
        source : `su.core.sources.Source`
            The source to load

        Returns
        -------
        rdts : dict
            A dictionary of region-dispersion tables (RDTs), where the
            keys are the region IDs and the values are the RDTs
        """
        rdts = {}
        for regid, region in enumerate(source):
            rdt = RDT(source, regid)

            for x, y in region.pixels():
                pdt = self.load_pdt(x, y)
                if len(pdt) == 0:
                    LOGGER.debug("MAJOR ERROR IN READING REGION")
                rdt.append(pdt)
            rdt.decimate()

            for k, v in pdt.attrs.items():
                if k not in ('x', 'y'):
                    rdt.attrs[k] = v

            rdts[regid] = rdt
        return rdts

    def load_pdts(self, source, flatten=False):
        """
        Method to load a bunch of pixel-dispersion tables (PDTs) from this
        `PDTFile`object-profile table (OPT) from this `PDTFile`
        given a source

        Parameters
        ----------
        source : `su.core.sources.Source`
            The source to load

        flatten : bool, optional
            Flag to flatten the collection of the PDTs.  Default is False


        Returns
        -------
        pdts : dict
            A dictionary of the PDTs.
            If flatten==True, then pdts will be a single dict, where the keys
               are just the (x,y) tuples.
            If flatten==False, then the pdts will be a nested dict, where
               the inner keys are for the Region ID and the outer are for
               the (x,y) tuples
        """

        pdts = {}
        if flatten:
            for x, y in source.pixels():
                xyd = source.image_coordinates(x, y, dtype=int)
                pdts[xyd] = self.load_pdt(*xyd)
        else:
            for regid, region in source.items():
                pdts[regid] = {}
                for x, y in region.pixels():
                    xyd = source.image_coordinates(x, y, dtype=int)
                    pdts[regid][xyd] = self.load_pdt(*xyd)

        return pdts

    def load_pdt(self, x, y):
        """
        Method to load the a pixel-dispersion table (PDT) from this `PDTFile`

        Parameters
        ----------
        x : int
           The x-coordinate from the direct image to load

        y : int
           The y-coordinate from the direct image to load

        Returns
        -------
        pdt : `su.core.tables.PDT`
           The output PDT
        """

        pdt = PDT(x, y)

        # check if this order is present in the PDTFile
        if pdt.name in self.h5order:
            # load the order
            hd = self.h5order[pdt.name]

            # load and copy over any attributes (think header keywords)
            for attr in hd.attrs:
                pdt.attrs[attr] = attributes.load(hd, attr)

            # copy the data into the columns of the output PDT
            data = hd[()]
            for column in pdt.COLUMNS:
                pdt[column] = data[column]

        return pdt
