from ...utilities import headers


class GroupCollection(list):
    """
    Class to collect the groups, which are sets of sources whose spectral
    traces overlap in a given set of WFSS images

    inherits from `list`

    Parameters
    ----------
    order : list or None
       The spectral orders to tabulate.  If `None`, then will do all orders
       present in the configuration file.  Default is None

    minarea : float, optional
        The minimum fractional area to test for overlaps, see also
        `su.core.modules.group.Group`  Default is 0.1


    """

    def __init__(self, orders=None, minarea=0.1):
        self.minarea = minarea
        self.orders = orders

    def __str__(self):
        if self:
            s = "Extraction Groups:\n" + '\n'.join('  '+str(g) for g in self)
        else:
            s = 'Extraction Groups (None)'
        return s

    def __call__(self, sources):
        if len(self) > 0:
            for group in self:
                yield {sources[g].segid: sources[g] for g in group}
        else:
            yield sources

    def append(self, grp):
        """
        Method to include a new group

        Parameters
        ----------
        grp : list or scalar
            The segmentation IDs for this group
        """

        super().append(tuple(grp))

    @classmethod
    def empty(cls):
        """
        Classmethod to create an empty `GroupCollection`.
        """
        return cls(minarea=None, orders=None)

    def find_group(self, segid):
        """
        Method to find the group ID for a given segmentation ID

        Parameters
        ----------
        segid : int
            The segmentation ID to search for

        Returns
        -------
        grpid : int
            The group ID that contains the input segmentation ID.  If the
            segmentation ID is not found, then returns None
        """

        for group, segids in enumerate(self):
            if segid in segids:
                return group

    def update_header(self, hdr):
        """
        Method to update a fits header

        Parameters
        ----------
        hdr : `astropy.io.fits.Header`
           The fits header to update
        """

        m = str(self.minarea) if self.minarea else 'N/A'
        o = ','.join(self.orders) if self.orders else "N/A"

        ngroup = max(len(self), 1)
        # grouped=ngroup>1
        # hdr['grouped']=(grouped,'Was this grouped?')
        hdr['groups'] = (ngroup, 'Number of groups')
        hdr['grparea'] = (m, 'Minimum fractional area for polygon intersection')
        hdr['grpord'] = (o, 'Orders checked for grouping')
        headers.add_stanza(hdr, 'Group Parameters', before='groups')
