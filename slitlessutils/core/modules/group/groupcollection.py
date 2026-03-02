import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nwx

from ....logger import LOGGER
from ...utilities import headers


class GroupCollection:
    """
    Use NetworkX to create graphs of spectral traces to assess
    contamination.  Objects are represented by nodes in the graph
    while contamination is represented by an edge, whose weight
    is specified when defining a contaminant

    Inputs
    ------
    threshold : float
        Threshold to flag contamination (not used, just for safe keeping)

    orders : str
        The orders used for contamination (not used, just for safe keeping)

    """

    def __init__(self, threshold=0.0, orders=''):

        # a basic graph to do the juggling
        self.graph = nwx.Graph()

        self.threshold = threshold
        self.orders = orders

    def add_object(self, segid, mag):
        """
        Add a spectral trace to the graph.

        Inputs
        ------
        segid : int, str, or tuple of int/str
            The Object ID to consider.

        mag : float
            The magnitude.  This is used to set color of node in
            graph.

        """
        # this is a basic node, but will add more attributes
        self.graph.add_node(segid, magnitude=mag)

    def add_contamination(self, segid0, segid1, weight):
        """
        Add a contamination entry

        Inputs
        ------
        segid0 : int, str or tuple of int/str
           The first object ID

        segid1 : int, str or tuple of int/str
           The second object ID

        weight : float
           The weight of the contamination.

        """
        # don't permit monitoring of self-contamination
        if segid0 == segid1:
            return

        if self.graph.has_edge(segid0, segid1):
            # if the edge exists, then update the weight
            self.graph[segid0][segid1]['weight'] += weight

        else:
            # if the edge doesn't exist, then make a new one
            self.graph.add_edge(segid0, segid1, weight=weight)

    def update_header(self, hdr, grouped=True):
        """
        Update a header for safe keeping

        Inputs
        ------
        hdr : `astropy.io.fits.Header`
           The header to update

        grouped : bool
           Flag if it was grouped (default is True)

        """

        hdr.set('GROUPED', value=grouped, comment='grouping flag')
        hdr.set('GRPTHRSH', value=self.threshold,
                comment='threshold for contamination')
        hdr.set('GRPORDER', value=self.orders,
                comment='orders for contamination')
        headers.add_stanza(hdr, 'Grouping Properties', before='GROUPED')

    def plot(self, filename=None, mag0=20, mag1=26, minobj=2):
        """
        Make a quick plot for diagnostics

        Inputs
        ------
        filename : str, optional
            A file name for output graphic.  Default is None

        mag0 : float, optional
            The minimum magnitude to assign to color mapping.
            Default is 20

        mag1 : float, optional
            The maximum magnitude to assign to color mapping.
            Default is 26

        minobj : int, optional
            The minimum number of nodes to show the subgraph.
            Default is 2

        """

        LOGGER.info(f"Making group graphic: {filename}.")

        # make copy for safekeeping
        graph = self.graph.copy()

        # remove distonnected subgraphs if they contain fewer than a
        # minimum number of nodes
        if minobj > 1:

            # collect the segids to remove
            badids = []
            for comp in nwx.connected_components(graph):
                if len(comp) < minobj:
                    badids.extend(comp)

            # remove them
            graph.remove_nodes_from(badids)

        # some plotting positions
        b, r = 0.03, 0.95
        w = (0.9, 0.02)
        l = (0.01, w[0] + 0.01)

        # make the plot window
        fig = plt.figure(figsize=(10, 8))
        gax = fig.add_axes((l[0], b, w[0], r))
        cax = fig.add_axes((l[1], b, w[1], r))

        # get positions for a given layout
        pos = nwx.spring_layout(graph)

        # Draw nodes with colormap
        mag = [a['magnitude'] for a in graph.nodes.values() if a]

        # draw and label the nodes
        nodes = nwx.draw_networkx_nodes(graph, pos, vmin=mag0,
                                        vmax=mag1, node_color=mag, cmap=plt.cm.coolwarm,
                                        node_size=200, ax=gax)
        nwx.draw_networkx_labels(graph, pos, font_size=10, ax=gax)

        # draw the edges.  Will scale the weights by some value and
        # introduce a minimum (don't want to have a weight that is
        # too small), but this scalar likely needs adjusting
        # based on how "weight" changes
        weight = [min(5 * graph[u][v]['weight'], 1) for u, v in graph.edges()]
        weight = [5 * graph[u][v]['weight'] for u, v in graph.edges()]
        nwx.draw_networkx_edges(graph, pos, width=weight,
                                ax=gax, edge_color='black')

        # make a legend
        weights = [1, 2, 3, 4, 5]
        lines = [Line2D([0], [0], color='black', lw=w) for w in weights]
        labels = [str(w / 5.) for w in weights]
        gax.legend(lines, labels, title='fractional area', loc='upper left')

        # dont draw the axes, since they have no meaning here
        gax.axis('off')

        # Add colorbar
        cbar = fig.colorbar(nodes, cax=cax)
        cbar.set_label("magnitude (mag)")

        # save to a file
        if isinstance(filename, str):
            plt.savefig(filename)
        else:
            plt.show()

    def groups(self, sources):
        """
        Iterator over the disconnected groups

        Inputs
        ------
        sources : `SourceCollection`
           The sources to sift through

        Returns
        -------
        grouped : dict
           The segids
        """
        n = len(self.graph)

        if n == 0:
            yield sources

        else:
            comps = nwx.connected_components(self.graph)

            for segids in sorted(comps, key=len, reverse=True):
                dct = {segid: sources[segid] for segid in segids}
                yield dct


if __name__ == '__main__':
    import numpy as np

    N = 100
    mags = np.random.uniform(low=20, high=25., size=N)

    groups = GroupCollection()

    for i, m in enumerate(mags, start=1):
        groups.add_object(i, m)

    for i in range(N):

        r = np.random.uniform()
        if r > 0.7:
            j = np.random.randint(low=1, high=N)
            groups.add_contamination(i, j, np.random.uniform())

            if np.random.uniform() > 0.2:
                groups.add_contamination(j, i, np.random.uniform())

    groups.add_object((1, 2), 23.)

    groups.plot(filename='t.svg')
