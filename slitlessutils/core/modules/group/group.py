import numpy as np
from shapely import geometry

from .groupcollection import GroupCollection
from ..module import Module
from ...tables import PDTFile


class Group(Module):
    """
    Module to group sources based on the WFSS imaging

    inherits from `su.core.modules.Module`

    Parameters
    ----------
    orders : list or None, optional
       The spectral orders to tabulate.  If `None`, then will do all orders
       present in the configuration file.  Default is None

    minarea : float, optional
       The fractional area in common to be consider an *overlap*.  Default
       is 0.1 (ie. >10% of the area of one object is common).

    Notes
    -----
    A group is defined here as the complete collection of sources whose
    spectral order will overlap in a collection of WFSS imaging.

    """

    DESCRIPTION = "Grouping WFSS"

    def __init__(self, minarea=0.1, orders=None, **kwargs):
        Module.__init__(self, self.group_wfss, postfunc=self.group_return,
                        **kwargs)

        self.minarea = minarea
        self.orders = orders

    def group_wfss(self, data, sources, **kwargs):
        """
        Method to group the WFSS images


        Parameters
        ----------
        data : `WFSSCollection`
           The WFSS data to group

        sources : `SourceCollection`
           The sources to group

        kwargs : dict, optional
           Dictionary of optional keywords.  Not currently used.

        Returns
        -------
        segids : list
           A list of sets, where each element of the list is a given
           group.  Therefore, this list has 1<=len(segids)<=len(sources).


        Notes
        -----
        This is the main routine, but is unlikely to be directly called

        """

        # conf,wfss=data
        with PDTFile(data, path=self.path, mode='r') as h5:
            groups = []
            for detname, detconf in data.items():
                h5.load_detector(detname)

                # the collection of Shapely Polygons and segids
                polys = []
                ids = []
                for i, (segid, source) in enumerate(sources.items()):
                    poly = []

                    # read each order as a Shapely Polygon
                    for ordname in self.orders:
                        h5.load_order(ordname)
                        odt = h5.load_odt(source)
                        poly.append(odt.shapelyPolygon())

                    # glue the individual orders together as a MultiPolygon
                    poly = geometry.MultiPolygon(poly)

                    # aggregate all the MultiPolygons and orders
                    polys.append(poly)
                    ids.append([source.segid])

                # package the ids and polygons
                data = list(zip(ids, polys))

                # group all of the sources' polygons using Shapely math
                grouped = self.group_polygons(data)

                # collect the grouped polygons
                groups.append(grouped)

        # at this point, we no longer need the shapely polygons. and instead
        # will just work with SEGIDs as sets, and use set math (is faster)
        segids = list(list(zip(*groups[0]))[0])
        segids = [set(i) for i in segids]
        segids = self.group_ids(segids)

        return segids

    def group_polygons(self, data):
        """
        Helper function to group sources based on shapely math

        Parameters
        ----------
        data : list
           A list of 2-tuples to group.  The elements of the tuple are given
           as the segmentation ID and a `shapely.geometry.Polygon` object
           respectively.

        Returns
        -------
        new : list
           A list of 2-tuples that have been grouped.  The output tuples
           have the same meaning as the input

        """
        nnew = ndata = len(data)

        # process until there is nothing left to process
        while nnew:
            groups = []

            while data:

                # grab the first polygon and ID
                thisid, thispoly = data.pop(0)

                # iterate over all other polygons and IDs
                for i, (testid, testpoly) in enumerate(data):
                    inter = thispoly.intersection(testpoly)

                    # is the intersection a positive region?
                    if (inter.area > self.minarea*testpoly.area) and \
                       (inter.area > self.minarea*thispoly.area):

                        # remove the test polygon, glue the test polygon and
                        # segIDs to the primary polygon and list
                        data.pop(i)
                        thispoly = thispoly.union(testpoly)
                        thisid.extend(testid)

                # collect this result
                groups.append((thisid, thispoly))

            # iterate until we've cycled thru everything
            N = len(data)
            data = groups
            nnew = ndata-N
            ndata = N
        return data

    def group_ids(self, data):
        """
        Helper method to group segmentation IDs using set-logic

        Parameters
        ----------
        data : list
            A list of sets.  Each element of the list represents a putative
            group, while each element of any set represents the segmentation
            ID of a source that has been grouped.

        Returns
        -------
        new : list
            A list of sets, with the same interpretation as the input list.
        """

        nnew = ndata = len(data)

        while nnew:
            new = []
            while data:
                this = data.pop(0)
                for i, test in enumerate(data):
                    if this.intersection(test):
                        this = this.union(test)
                        data.pop(i)
                new.append(this)
            data = new
            n = len(data)
            nnew = ndata-n
            ndata = n
        return data

    def group_return(self, ids, data, sources, **kwargs):
        """
        Helper method to package the results into a `GroupCollection`

        Parameters
        ----------
        ids : list
            The list of sets from `self.group_ids()`

        data : `WFSSCollection`
            The WFSS data that was used to group

        sources : `SourceCollection`
            The sources to group

        kwargs : dict, optional
            Not used, but is a place holder for optional parameters

        Returns
        -------
        out : `GroupCollection`
           A collection of groups

        """

        sets = []
        for i in ids:
            sets.extend(i)

        # group the IDs and sort by size
        groups = self.group_ids(sets)
        groups.sort(key=len, reverse=True)

        # sort out what is the output product going to look like
        out = GroupCollection(minarea=self.minarea, orders=self.orders)
        for grp in groups:
            out.append(grp)

        # out={k:list(v) for k,v in enumerate(groups)}
        # out=[list(g) for g in groups]

        return out
