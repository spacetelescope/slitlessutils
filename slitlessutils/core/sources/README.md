# Sources
This is a collection of routines for working with *sources*.  The main data structure is a *source collection* that is built from a *segmentation map*.  The *source collection* will operate like a python `dict`, and a *source* will operate like a python `list`.

An example source collection with 3 sources and several spectral regions (7 in this example).

```
source collection
|
+-- Source #1
|   |
|   +-Spectral region 1.1
|   +-Spectral region 1.2
|
+-- Source #2
|   |
|   +-Spectral region 2.1
|
+-- Source #3
    |
    +-Spectral region 3.1
    +-Spectral region 3.2
    +-Spectral region 3.3
    +-Spectral region 3.4

```
Where each spectral region has a unique spectrum, that is described a structured `np.ndarray` (see the `photometry` sub-package).

## Glossary

**segmentation map:** a 2d image that describes where the *sources* or *spectral regions* are on the sky.  In the simplest case, this can directly taken as the output of [Source Extractor](https://sextractor.readthedocs.io/en/latest/Introduction.html), but more complex options exist.  See the readthedocs for more details.

**source:** an astrophysical object or region that has some coherent morphological structure, but can be decomposed into many *spectral regions*.  This will act as a ```list```, where the elements are spectral regions.

**source collection:** a compendium of sources, that is generally instantiated from a *segmentation map*.  This will act as as a ```dict```, whose keys are the segmentation IDs and values are a ```Source()```.

**spectral region:** a subset of a source, where the spectrum is assumed to be constant.  But this object records the region (in the form of the direct-image pixels) correspond to this region.  A source is composed of one or more spectral regions.
