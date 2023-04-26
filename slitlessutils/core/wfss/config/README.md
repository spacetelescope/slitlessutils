# WFSS Config
A series of routines for configuring the WFSS data.  In general, many of these are ```dataclasses``` to configure various subcomponents.


## Main Configuration Objects


### **wfssconfig.py:** 

This configures the spectrscopic setup.  This is equivalent to ```grismconf``` or ```aXeConfig```, and will act as a ```dict``` where the keys are the spectral-order names and the values are from the ```order.py```.  The spectroscopic orders contain the primary trace and dispersion objects, 

- ```dispx``` and ```dispy``` represent the $x$ and $y$ coordinates of the spectral trace (respectively), and are enocded as ```StandardPolynomial``` objects in the ```parametricpolynomial.py``` file.  A "standard polynomial" is given as:
$$x(t;x_0,y_0) = a_0(x_0,y_0) + a_1(x_0,y_0) t + a_2(x_0,y_0)t^2 ...$$
$$y(t;x_0,y_0) = b_0(x_0,y_0) + b_1(x_0,y_0) t + b_2(x_0,y_0)t^2 ...$$
where the functions $\{a(x_0,y_0)\}$ and $\{b(x_0,y_0)\}$ are given as ```SpatialPolynomial```s:
$$a_i(x_0,y_0) = a_{0,0} + a_{0,1} x_0 + a_{0,2} y_0 + a_{0,3} x_0^2 + a_{0,4}x_0y_0 + a_{0,5}y_0^2 + ...$$
where the coefficients $\{a\}$ are given as calibrations (and similarly for $\{b\}$).  These are given as an array, whose length is a [triangular number](https://en.wikipedia.org/wiki/Triangular_number) and encoded as [Cantor pairs](https://en.wikipedia.org/wiki/Pairing_function).  *nota bene:* this representation is the same as the SIP coefficients.  

-  ```displ``` represents the wavelength along the trace.  For the grism spectral elements (which generally start with either "G" or "GR"), this will be encoded as a ```StandardPolynomial```.  However, the prism spectral elements (which generally start with "P" or "PR") will be encoded as ```ReciprocalPolynomial()``` objects in the ```parametricpolynomial.py``` file.  These are given as:
$$\lambda(t;x_0,y_0) = c_0(x_0,y_0) + \frac{c_1(x_0,y_0)}{(t-t^*(x_0,y_0))}+ \frac{c_2(x_0,y_0)}{(t-t^*(x_0,y_0))^2}+...$$
where the $\{c\}$ and $t^*$ functions are given as ```SpatialPolynomials```.

In any case, the process of dispersing a region results in a 3-step procedure:

1. Take direct image position $(x,y)$ and use WCS transforms to convert to an "undispersed" position $(x_0,y_0)$ for a given ```WFSSDetector```.
2. Specify a wavelength, and compute the $t$ value from the inverse of ```displ```:  
$$t = DISPL^{-1}(\lambda;x_0,y_0)$$
3. These parameters are passed into the trace functions to give the grating-specific position on a ```WFSSDetector```:
$$x_g=x(t;x_0,y_0)$$
$$y_g=y(t;x_0,y_0)$$

### **instrumentconfig.py**

This configures the instrumental setup, which describes many things.  

- The instrument: ```InstrumentConfig```.  This acts as ```dict```, where the keys and values are the detector names and ```DetectorConfig``` objects.  But, unlike a ```dict```, it further contains information like the properties of the disperser and the relative positions of the detectors on the focal plane (also called the ```SIAF``` data).  
- The ```DetectorConfig``` governs a given detector, and contains the pixelation properties (ie. CRPIX, pixel scales, image sizes, bad-pixel values, and additional header information), the file format properties (e.g. the fits MEF information), and some science-based data (such as the noise model).  
- The ```SIAF``` dataclass provides the relative positions of the detectors in the focal plane.  ***THIS IS ONLY USED FOR SIMULATIONS.***
- The ```Extension``` dataclass describes a MEF extension.

### Miscellaneous Configuration Objects
- **disperser.py** this governs the grating and its dispersive properties (eg resolution, range, etc.)
- **flatfield.py** this describes a WFSS flat field, which can be either:
  -  ```UnityFlatField``` assume the flat is all 1.0 (effectively no flat fielding)
  -  ```ImageFlatField``` this is effectively a "gray" flat that is generally taken from an image flat.
  -  ```PolynomialFlatField``` this is true, wavelength-dependent flat field described as a polynomial in wavelength:
   $$F(\lambda,x_g,y_g) = F_0(x_g,y_g) + F_1(x_g,y_g)\left(\frac{\lambda-\lambda_0}{\lambda_1-\lambda_0}\right) + F_2(x_g,y_g)\left(\frac{\lambda-\lambda_0}{\lambda_1-\lambda_0}\right)^2...$$
where $\lambda_0$ and $\lambda_1$ are the range of the grating (and specified in the calibration), and the images $\{F\}$ are given in a MEF file from calibration.
- **pom.py** This describes the "pick-off mirror" and the effective size of the field-of-regard (ie. the region of the sky that surrounds a ```WFSSDetector``` that could potentially disperse light onto the detector), which can be either:
  - ```UnityPOM``` assume the POM is all 1.0 (effectively no POM restrictions)
  - ```RangePOM``` a POM that is 1.0 within a rectangular range and 0.0 outside of it.  This is typical for the HST instruments.
  - ```ImagePOM``` a POM that can be a floating point value, to restrict a particular *fraction* of the incident light.  This is common for the JWST instruments
- **sensitivity.py** inherits from the ```photometry``` subpackage and holds the information for the sensitivity curve