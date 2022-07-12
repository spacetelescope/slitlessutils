# Astrometry Module Notes

These are very old classes that need a bit of spit-and-polish, and they may be totally redundant with other tools.

1) ```astroimage``` facilitates a few other operations, including cutouts, smoothing, and rebinning.  Some of this might replicated in some combination of astropy and scikit-image, but this works for the moment.

2) ```wcs``` extends the ```astropy.wcs.WCS``` object with a few more calculations and convience functions. 