# vectorialpsf
Python bindings for a vectorial model of the point spread function implemented in C++. 

Original C code Copyright &copy; 2006-2013, [Francois Aguet](https://github.com/francois-a), distributed under GPL-3.0 license. 
Python bindings by Talley Lambert, 2019.

The model is described in Auget et al 2009<sup>1</sup>.  For more information and implementation details, see Francois' Thesis<sup>2</sup>.

<sup>1</sup> [F. Aguet et al., (2009) Opt. Express 17(8), pp. 6829-6848](https://doi.org/10.1364/OE.17.006829)

<sup>2</sup> [F. Aguet. (2009) Super-Resolution Fluorescence Microscopy Based on Physical Models. Swiss Federal Institute of Technology Lausanne, EPFL Thesis no. 4418 (2009), 209 p., May 28, 2009.](http://bigwww.epfl.ch/publications/aguet0903.html)

### see also:

For a different (faster) scalar-based Gibson–Lanni PSF model, see the [MicroscPSF](https://github.com/MicroscPSF) project, based on [Li et al (2017)](https://doi.org/10.1364/JOSAA.34.001029) which has been implemented in Python, MATLAB, and ImageJ/Java

## Usage

There are two main functions in `vectorialpsf`: `make_psf` and `make_centered_psf`.  The main difference is that make_psf accepts a vector of Z positions `zvec` (relative to coverslip) at which PSF is calculated, and an argument `zdepth`, specifying the position of the point source relative to the coverslip.  As such, the point source may or may not actually be in the center of the rendered volume.  `make_centered_psf`, by contrast, does not accecpt `zvec`, but rather accepts `nz` (the number of z planes) and `range` (the range in microns of the full zplane), and always generates an output volume in which the point source is positioned in the middle of the Z range, with planes equidistant from each other and a Z step size of `range / (nz - 1)`.

*Note, all output dimensions (`nx` and `nz`) should be odd.*

```python
import vectorialpsf as vpsf
import numpy as np

psf = vpsf.make_centered_psf(nz=63, nx=127)
psf.shape  # (63, 127, 127)

zvec = np.linspace(-3e-6, 3e-6, 31)
psf = vpsf.make_psf(zvec, nx=127)
psf.shape  # (31, 127, 127)
```

Both functions accept the following parameters.  In general, units should be provided in microns.

```
nx (int):       XY size of output PSF in pixels, must be odd.
zdepth (float): depth of point source relative to coverslip (in microns)
ti0 (float):    (optional, default: 1.515)  working distance of the objective (microns)
ni0 (float):    (optional, default: 1.518) immersion medium refractive index, design value
ni (float):     (optional, default: 1.518) immersion medium refractive index, experimental value
tg0 (float):    (optional, default: 170) coverslip thickness, design value (microns)
tg (float):     (optional, default: 170) coverslip thickness, experimental value (microns)
ng0 (float):    (optional, default: 1.515) coverslip refractive index, design value
ng (float):     (optional, default: 1.515) coverslip refractive index, experimental value
ns (float):     (optional, default: 1.47) sample refractive index
lambda (float): (optional, default: 0.55) emission wavelength (microns)
mag (float):    (optional, default: 100) magnification
na (float):     (optional, default: 1.45) numerical aperture
pixel (float):  (optional, default: 6.5) physical size (width) of the camera pixels (microns)
sf (int):       (optional, default: 3) oversampling factor to approximate pixel integration
mode (int):     (optional, default: 1) if 0, returns oversampled PSF
```