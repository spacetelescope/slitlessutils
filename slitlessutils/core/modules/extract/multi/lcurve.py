import numpy as np
import datetime
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
import matplotlib.colors as mc
from matplotlib.backends.backend_pdf import PdfPages

from .menger import menger


class LCurve:
    """
    Class to contain data for an L-curve analysis.

    Parameters
    ----------
    norm : float, optional
        Normalization constant for the ||x||^2 term.  In Ryan+ 2018, this
        is taken as the Frobenius norm of the linear operator :math:`A`.
        Default is 1.0 (ie. no normalization)

    Notes
    -----
    1) An L-curve

    2) Internally, the data is stored in an `astropy.table.Table`

    """

    # variables for writing L-Curve data to fits tables
    NAMES = ('iter', 'logchi2', 'lognorm', 'logdamp', 'curvature')
    DTYPES = (int, float, float, float, float)
    FORMATS = ('4d', '+.4f', '+.8e', '+.8e', '+.8e')

    def __init__(self, norm=1.):
        self.norm = norm
        self.clear()

    def __len__(self):
        return len(self.data)

    def __str__(self):
        return str(self.data)

    def clear(self):
        """
        Method to clear this L-curve object from all its data
        """

        self.it = 0
        self.data = Table(names=self.NAMES, dtype=self.DTYPES)

    def append(self, chi2, xnorm, logl):
        """
        Method to append a new row onto the table

        Parameters
        ----------
        chi2 : float
           The chi2 of a given fit

        xnorm : float
           The R2 norm of x (ie. the sum of the squares)

        logl : float
           The log10 of the damping parameter
        """

        self.it += 1
        self.data.add_row((self.it, np.log10(chi2), np.log10(xnorm), logl, np.nan))

    def compute_curvature(self):
        """
        Method to compute the local curvature at every point.

        Notes
        -----
        1) Curvature is defined as the Menger curvature (see menger.py
           and/or https://en.wikipedia.org/wiki/Menger_curvature)
        2) Because the Menger curvature requires a point before and a
           point after, the end points of the curvature array will be
           `np.nan`.

        """

        self.data.sort('logchi2')
        for j in range(1, len(self.data)-1):
            i = j-1
            k = j+1

            pi = (self.data['logchi2'][i], self.data['lognorm'][i])
            pj = (self.data['logchi2'][j], self.data['lognorm'][j])
            pk = (self.data['logchi2'][k], self.data['lognorm'][k])
            self.data['curvature'][j] = menger(pi, pj, pk)

    def as_HDU(self, grpid=0, **kwargs):
        """
        Method to package this curvature table as a header-data unit (HDU)

        Parameters
        ----------
        grpid : int, optional
           The group ID for this L-curve data.  Default is 0.

        kwargs : dict, optional
           Additionoal keyword/value pairs added as header cards
        """

        self.compute_curvature()

        hdu = fits.BinTableHDU(data=self.data)
        hdu.header['EXTNAME'] = ("LCURVE", 'extension name')
        hdu.header['EXTVER'] = (grpid, 'extension version number')
        hdu.header['GRPID'] = (grpid, 'Group ID')
        hdu.header['FROBNORM'] = (self.norm, 'Frobenius norm')

        for k, v in kwargs.items():
            hdu.header[k] = v

        return hdu

    def plot(self, arg, **kwargs):
        """
        Helper method to make an L-curve plot

        Parameters
        ----------
        arg : str or `matplotlib.pylot.PdfPages`
           if `str`, then will be the name of the pdffile
           if `PdfPages`, then will be a PDF file to append a new plot to

        kwargs : dict, optional
           Additional arguments set to `self.pdfplot()`
        """

        if isinstance(arg, str):
            with PdfPages(arg) as pdf:
                self.pdfplot(pdf, **kwargs)
        elif isinstance(arg, PdfPages):
            self.pdfplot(arg, **kwargs)
        else:

            LOGGER.error(f"Invalid argument type in plot {arg}")

    def pdfplot(self, pdf, colormap='Spectral', lighten=0.5, grpid=None):
        """
        Method to plot the L-curve

        Parameters
        ----------
        pdf : `matplotlib.pyplot.PdfPages`
           A pdf plotting device

        colormap : str, optional
           The color map to plot with.  Default is 'Spectral'

        lighten : float, optional
           A parameter to lighten the colormap, must be 0<=lighten<=1.
           Default is 0.5.

        grpid : int, optional
           The group ID --- used only to add to plot.  Default is None

        """

        # get the colormap and lighten it
        cmap = plt.cm.get_cmap(colormap)
        cmap = self.remap_cmap(lambda x: lighten*(1.+x), cmap)

        # compute curvatures
        self.compute_curvature()

        # some plot windows
        fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [0.8, 1.]})

        # show the group ID
        if grpid is not None:
            ax1.set_title(f'Group ID: {grpid}')

        lmin = np.amin(self.data['logdamp'])
        lmax = np.amax(self.data['logdamp'])

        # create dummy variables for the plotting
        logchi2 = self.data['logchi2']
        lognorm = self.data['lognorm']

        # scale the norm if need be
        rescale = self.norm is not None
        if rescale:
            lognorm += np.log10(self.norm)

        # plot the data in black and data points in color
        #line = ax1.plot(self.data['logchi2'], self.data['lognorm'], '-k',
        line = ax1.plot(logchi2, lognorm, '-k', linewidth=1, zorder=1)

        #scat = ax1.scatter(self.data['logchi2'], self.data['lognorm'], zorder=2,
        scat = ax1.scatter(logchi2, lognorm, zorder=2,
                           c=self.data['logdamp'], s=40, cmap=cmap, edgecolors='k',
                           marker='o', vmin=lmin, vmax=lmax)

        # put the frobenius norm in the plot
        #if self.norm is not None:
        if rescale:
            text = ax1.text(0.67, 0.88,
                            r'$\log\ ||W||_F=${0:+.3f}'.format(self.norm),
                            horizontalalignment='left', transform=ax1.transAxes,
                            bbox=dict(facecolor='white', edgecolor='white'))
            ylabel = r'$\log\ ||W||_F^2\,||f||^2$'
        else:
            ylabel = r'$\log\ ||f||^2$'

        # plot the labels
        #ax1.set(xlabel=r'$\log\ ||Ax-b||^2$',
        #        ylabel=r'$\log\ ||x||^2$')
        ax1.set(xlabel=r'$\log\ ||Wf-{\cal I}||^2$',
                ylabel=ylabel)
        ax1.set_axisbelow(True)

        # put on a grid
        ax1.grid(True)

        # plot the curvature

        line = ax2.plot(self.data['logdamp'], self.data['curvature'], '-k',
                        linewidth=1, zorder=1)

        scat = ax2.scatter(self.data['logdamp'], self.data['curvature'],
                           zorder=3, s=40, cmap=cmap, edgecolors='k', marker='o',
                           c=self.data['logdamp'], vmin=lmin, vmax=lmax)

        # find and indicate the maximum
        if len(self.data) > 1:
            g = np.where(self.data['curvature'] == np.nanmax(self.data['curvature']))[0][0]
        else:
            g = 0
            ax2.axvline(x=self.data['logdamp'][g], color='k', linestyle=':', zorder=2)

        # set the xrange
        if lmax == lmin:
            if np.isfinite(lmax):
                ax2.set_xlim([lmin-0.5, lmax+0.5])
            else:
                ax2.set_xlim([-5, 0])
        else:
            ax2.set_xlim([lmin, lmax])

        # label the axes
        ax2.set(ylabel='curvature')
        ax2.set_xticklabels([])
        ax2.set_axisbelow(True)

        # make a colobar
        cbar = plt.colorbar(scat, ax=ax2, orientation='horizontal', pad=0.0,
                            aspect=60, shrink=1)
        cbar.ax.set_xlabel(r"$\log\ \ell$")

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()

    @staticmethod
    def remap_cmap(function, cmap):
        """
        Staticmethod to remap a color map.

        Parameters
        ----------
        function : callable
            The function to apply to the colormap

        cmap : `matplotlib.colors.ColorMap`
            The colormap to remap

        Notes
        -----
        1) Will breka any discontinuous points in a colormap
        2) operates on RGB vectors
        """
        cdict = cmap._segmentdata
        step_dict = {}
        # Firt get the list of points where the segments start or end
        for key in ('red', 'green', 'blue'):
            step_dict[key] = list(map(lambda x: x[0], cdict[key]))
        step_list = sum(step_dict.values(), [])
        step_list = np.array(list(set(step_list)))
        # Then compute the LUT, and apply the function to the LUT
        def reduced_cmap(step): return np.array(cmap(step)[0:3])
        old_LUT = np.array(list(map(reduced_cmap, step_list)))
        new_LUT = np.array(list(map(function, old_LUT)))
        # Now try to make a minimal segment definition of the new LUT
        cdict = {}
        for i, key in enumerate(['red', 'green', 'blue']):
            this_cdict = {}
            for j, step in enumerate(step_list):
                if step in step_dict[key]:
                    this_cdict[step] = new_LUT[j, i]
                elif new_LUT[j, i] != old_LUT[j, i]:
                    this_cdict[step] = new_LUT[j, i]
            colorvector = list(map(lambda x: x + (x[1], ), this_cdict.items()))
            colorvector.sort()
            cdict[key] = colorvector

        return mc.LinearSegmentedColormap('colormap', cdict, 1024)

    def write(self, filename, **kwargs):
        """
        Method to write the Lcurve data to a file

        Parameters
        ----------
        filename : str
           the filename to write.  The type of file is taken from the
           extension in this file.  Extensions of: txt, ascii, dat, lcv,
           and/or lcurve are interpeted as ascii files, see
           `self.write_ascii()` for more information.  Extensions of
           fit or fits are fits files, see `self.write_fits()`

        kwargs : dict, optional
           Additional keywords passed to individual writers.

        """

        ext = os.path.splitext(filename)[-1]
        if ext in ('txt', 'ascii', 'dat', 'lcv', 'lcurve'):
            self.write_ascii(filename, **kwargs)
        elif ext in ('fits', 'fit'):
            raise NotImplementedError("fits files not yet supported")

            # self.write_fits(filename,**kwargs)
        else:
            raise NotImplementedError(ext)

    def write_ascii(self, filename, delim=' '):
        """
        Method to write an ascii file

        Parameters
        ----------
        filename : str
            The filename to write

        delim : str, optional
            The delimiter for the file.  Default is ' '
        """

        now = datetime.datetime.now()
        date = now.strftime("%Y-%m-%d")

        x, y, l, s = self.values()
        c = self.compute_curvature()

        with open(filename, 'w') as fp:
            print(f'# File written by {__code__}: {date}')
            print(f'# FROBNORM = {self.norm}', file=fp)
            print('# 1: iteration', file=fp)
            print('# 2: log(ell)', file=fp)
            print('# 3: log(||Ax-b||^2)', file=fp)
            print('# 4: log(||x||^2)', file=fp)
            print('# 5: curvature', file=fp)
            for values in zip(s, l, x, y, c):
                out = [f'{{0:{f}}}'.format(v) for f, v in zip(self.FORMATS, values)]
                print(delim.join(out), file=fp)
