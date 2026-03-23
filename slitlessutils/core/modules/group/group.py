from collections import namedtuple

import networkx
import numpy as np
from shapely import Polygon

from ...tables import PDTFile
from ..module import Module
from .groupcollection import GroupCollection

Key = namedtuple('Key', ['filename', 'detname', 'ordname', 'segid'])


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

    def __init__(self, threshold=0.01, orders=None, **kwargs):

        Module.__init__(self, self.group_oneimage, postfunc=self.group_images,
                        **kwargs)

        self.threshold = np.clip(threshold, 0, 1)
        self.orders = orders

    def group_oneimage(self, data, sources, **kwargs):

        group = GroupCollection(threshold=self.threshold,
                                orders=self.orders)

        # read all sources as Shapely polygons
        polys = []
        with PDTFile(data, path=self.path, mode='r') as h5:
            for detname, detconf in data.items():
                h5.load_detector(detname)

                for ordname in self.orders:
                    h5.load_order(ordname)

                    for segid, source in sources.items():
                        odt = h5.load_odt(source)
                        if odt:

                            # load the vertices
                            x0, x1, y0, y1 = odt.bounding_box()
                            vertices = ((x0, y0), (x0, y1), (x1, y1),
                                        (x1, y0), (x0, y0))

                            # store some data for this trace
                            key = Key(data.filename, detname, ordname, segid)
                            poly = Polygon(vertices)

                            # save the trace data
                            polys.append((key, poly))

        # add vertices for each source
        for segid, source in sources.items():
            group.add_object(segid, source.mag)

        # now test all n^2 operations.  But we can ignore the lower triangle
        n = len(polys)
        for i in range(n):
            keyi, polyi = polys[i]

            for j in range(i + 1, n):
                keyj, polyj = polys[j]

                inter = polyi.intersection(polyj)
                if inter:
                    union = polyi.union(polyj)

                    ratio = inter.area / union.area
                    if ratio > self.threshold:
                        group.add_contamination(keyi.segid, keyj.segid, ratio)

        return group

    def group_images(self, groups, data, sources, **kwargs):
        group = GroupCollection()
        for grp in groups:
            group.graph = networkx.compose(group.graph, grp.graph)

        return group
