# WFSS Data
This is collection of routines for working with WFSS data (simulated or observed). The main data structure is the *wfss collection* 

```
wfss collection
|
+-- WFSS #1
|   |
|   +-- detector #1
|   +-- detector #2
+-- WFSS #2
|   |
|   +-- detector #1
|   +-- detector #2
+-- WFSS #3
|   |
|   +-- detector #1
|   +-- detector #2
+-- WFSS #4
|   |
|   +-- detector #1
|   +-- detector #2
+-- WFSS #5
    |
    +-- detector #1
    +-- detector #2
```








## Glossary

**blocking filter:** an imaging filter inserted in the optical path along with the *spectral element*.  This is only relevant to JWST, as none of the HST instruments support blocking filters --- as such these will be listed as either ```None``` or ```''``` (ie the empty string).  By and large, blocking filters are only needed to restrict the wavelength range, which minimizes the contamination.  


**ObservedData:** a ```dataclass``` that contains the metadata for an observed dataset.

**SimulatedData:** a ```dataclass``` that contains the metadata for a simulated dataset.

**spectral element:** the grating that sets one term of the spectral resolution.  This is usually specified as either a string or a 2-tuple; the first element will be a string of the grating name and the second is the blocking filter.


**WFSS:** wide-field slitless spectroscopy, but here is a short-hand for a wide-field slitless spectroscopic image.  This will act as a ```dict```, whose keys are detector names and values are 

**WFSS collection:** a compendium of WFSS images.  In general, these can be quite large in memory usage, therefore, only the light-weight metadata is loaded (e.g. WCS, filter names, etc.).  This will act as a ```dict```, whose keys are the dataset name (ie. IPPPSSOOT) and the values are either an ```ObservedData``` or ```SimulatedData``` as relevant.  

**WFSS Detector:** a given instrument may be composed of multiple detectors (e.g. WFC3/UVIS has 2 CCDs).  This is a key component, as this is where the WFSS data are actually held --- but as with the ```WFSSCollection()```, only the metadata are actually stored.  There are methods to read the pixel values.  

**WFSS Instrument:** this is a nebulously defined thing, as there is no single definition that satisfactorily works for all telescopes and "instruments".  Therefore, here we will adopt a pragmatic viewpoint, where an "instrument" is a single-setup or usage.  The currently supported instruments are given in this 

| Telescope |  Instrument | Spectral elements | 
|-----------|-------------|-------------------|
| HST       | WFC3IR      | G102, G141        |
|           | WFC3UVIS    | G280              |
|           | ACSWFC      | G800L             |
|           | ACSSBC      | P110L, P130L      |
| JWST      | NIRISS      | GR150R<sup>&dagger;</sup>, GR150C<sup>&dagger;</sup> |

<sup>&dagger;</sup>The JWST instruments generally have a grating that has a very broad wavelength range, and users "select" a subset of this using a *blocking filter*.  
