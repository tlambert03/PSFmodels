try:
    import _psfmodels as library
except ImportError:
    raise ImportError("could not import precompiled library for psfmodels")
import numpy as _np
import warnings as _wrn


_DEFAULT_PARAMS = {
    "NA": 1.4,  # numerical aperture
    "ng0": 1.515,  # coverslip RI design value
    "ng": 1.515,  # coverslip RI experimental value
    "ni0": 1.515,  # immersion medium RI design value
    "ni": 1.515,  # immersion medium RI experimental value
    "ns": 1.47,  # specimen refractive index (RI)
    "ti0": 150,  # microns, working distance (immersion medium thickness) design value
    "tg": 170,  # microns, coverslip thickness experimental value
    "tg0": 170,  # microns, coverslip thickness design value
}

_VALID_ARGS = [
    "zv",
    "nx",
    "pz",
    "ti0",
    "ni0",
    "ni",
    "tg0",
    "tg",
    "ng0",
    "ng",
    "ns",
    "wvl",
    "NA",
    "dxy",
    "sf",
    "mode",
]


def normalize_params(mp):
    """Check and return valid microscope parameters dict, stripped of excess keys.

    Args:
        mp (dict): The microscope parameters dict

    Raises:
        ValueError: If one of the parameter values is invalid

    Returns:
        dict: Valid parameters for use in vectorial_ or scalar_psf
    """
    _mp = _DEFAULT_PARAMS.copy()
    if mp is not None:
        if isinstance(mp, dict):
            _mp.update(mp)
        else:
            raise ValueError("mp argument must be dict of microscope params")
    out = {}
    for key, val in _mp.items():
        if key in _VALID_ARGS:
            out[key] = val
        else:
            _wrn.warn(
                (
                    "parameter {} is not one of the recognized "
                    "keywords and is being ignored"
                ).format(key)
            )

    return {key: val for key, val in _mp.items() if key in _VALID_ARGS}


def _validate_args(zv, dxy, pz):
    """sanity checks for various arguments"""
    if isinstance(zv, (float, int)):
        zv = _np.array([zv])
    elif isinstance(zv, (list, tuple)):
        zv = _np.array(zv)
    elif not isinstance(zv, _np.ndarray):
        raise ValueError("zd must be a scalar, iterable, or numpy array")
    if not dxy > 0:
        raise ValueError("dxy must be greater than 0")
    if not pz >= 0:
        raise ValueError("pz should be >= 0")
    return zv


_infostring = """

    For more information and implementation details, see:
    F. Aguet et al., Opt. Express 17(8), pp. 6829-6848, 2009
    F. Aguet, Ph.D Thesis, Swiss Federal Institute of Technology, EPFL (2009)"""

_docstring = (
    _infostring
    + """

    Args:
        zv (float, list, np.ndarray): Float or vector of Z positions
            (relative to coverslip) at which PSF is calculated. Defaults to [0].
            len(zv) should be odd"""
)

_centerdocstring = (
    _infostring
    + """

    Args:
        nz (int): Z size of output PSF in pixels, must be odd.
        dz (float, optional): Z step size of PSF in sample space. Defaults to 0.05
    Kwargs:"""
)

_paramdocs = """
        nx (int, optional): XY size of output PSF in pixels, must be odd. Defaults to 31.
        dxy (float, optional): XY Pixel size in sample space (microns). Defaults to 0.05.
        pz (float, optional): Depth of point source relative to coverslip in microns.
                              Defaults to 0.
        wvl (float, optional): Emission wavelength in microns. Defaults to 0.6.
        params (dict, optional): Microscope parameters dict. See optional keys below.
        normalize (bool, optional): Normalize PSF peak to 1. Defaults to True.

        valid params (all floats unless stated, all distances in microns):
            NA:  Numerical Aperture. Defaults to 1.4
            ng0: Coverslip RI design value. Defaults to 1.515
            ng:  Coverslip RI experimental value. Defaults to 1.515
            ni0: Immersion medium RI design value. Defaults to 1.515
            ni:  Immersion medium RI experimental value. Defaults to 1.515
            ns:  Specimen refractive RI. Defaults to 1.47
            ti0: Working distance (immersion medium thickness) design value.
                 Defaults to 150
            tg:  Coverslip thickness experimental value. Defaults to 170.
            tg0: Coverslip thickness design value. Defaults to 170.
            sf (int):  oversampling factor to approximate pixel integration. Defaults to 3
            mode (int): if 0, returns oversampled PSF. Defaults to 1

    Returns:
        np.ndarray: The PSF volume.
 """


def vectorial_psf(zv=0, nx=31, dxy=0.05, pz=0.0, wvl=0.6, params=None, normalize=True):
    """Computes a vectorial model of the microscope point spread function."""
    zv = _validate_args(zv, dxy, pz)
    params = normalize_params(params)
    _psf = library.vectorial_psf(zv.copy(), int(nx), pz, wvl=wvl, dxy=dxy, **params)
    if normalize:
        _psf /= _np.max(_psf)
    return _psf


def scalar_psf(zv=0, nx=31, dxy=0.05, pz=0, wvl=0.6, params=None, normalize=True):
    """Computes the scalar PSF model described by Gibson and Lanni."""
    zv = _validate_args(zv, dxy, pz)
    params = normalize_params(params)
    _psf = library.scalar_psf(zv.copy(), int(nx), pz, wvl=wvl, dxy=dxy, **params)
    if normalize:
        _psf /= _np.max(_psf)
    return _psf


def _centered_zv(nz, dz, pz):
    lim = (nz - 1) * dz / 2
    return _np.linspace(-lim + pz, lim + pz, nz)


def vectorial_psf_centered(nz, dz=0.05, **kwargs):
    """Computes a vectorial model of the microscope point spread function.

    The point source is always in the center of the output volume."""
    zv = _centered_zv(nz, dz, kwargs.get("pz", 0))
    return vectorial_psf(zv, **kwargs)


def scalar_psf_centered(nz, dz=0.05, **kwargs):
    """Computes the scalar PSF model described by Gibson and Lanni.

    The point source is always in the center of the output volume."""
    zv = _centered_zv(nz, dz, kwargs.get("pz", 0))
    return scalar_psf(zv, **kwargs)


vectorial_psf.__doc__ += _docstring + _paramdocs
scalar_psf.__doc__ += _docstring + _paramdocs
vectorial_psf_centered.__doc__ += _centerdocstring + _paramdocs
scalar_psf_centered.__doc__ += _centerdocstring + _paramdocs


def vectorialXYZFocalScan(mp, dxy, xy_size, zv, normalize=True, pz=0.0, wvl=0.6, zd=0):
    """Computes a vectorial model of the microscope point spread function.

    This function is merely here as a convenience to mimic the MicroscPSF-Py API.
    """
    return vectorial_psf(zv, xy_size, dxy, pz, wvl, mp, normalize)


def scalarXYZFocalScan(mp, dxy, xy_size, zv, normalize=True, pz=0.0, wvl=0.6, zd=0):
    """Computes the scalar PSF model described by Gibson and Lanni.

    This function is merely here as a convenience to mimic the MicroscPSF-Py API.
    """
    return scalar_psf(zv, xy_size, dxy, pz, wvl, mp, normalize)


_otherapidoc = """
    Args:
        mp (dict): Microscope parameters dict. See optional keys below.
        dxy (float): XY Pixel size in sample space (microns).
        xy_size (int): XY size of output PSF in pixels, must be odd.
        zv (float, list, np.ndarray): Float or vector of Z positions
            (relative to coverslip) at which PSF is calculated. Defaults to [0].
            len(zv) should be odd.
        normalize (bool, optional): Normalize PSF peak to 1. Defaults to True.
        pz (float, optional): Depth of point source relative to coverslip in microns.
                              Defaults to 0.
        wvl (float, optional): Emission wavelength in microns. Defaults to 0.6.
        zd (int, optional): Unused in this module.  Just here for MicroscPSF-Py API.

        valid params (all floats unless stated, all distances in microns):
            NA:  Numerical Aperture. Defaults to 1.4
            ng0: Coverslip RI design value. Defaults to 1.515
            ng:  Coverslip RI experimental value. Defaults to 1.515
            ni0: Immersion medium RI design value. Defaults to 1.515
            ni:  Immersion medium RI experimental value. Defaults to 1.515
            ns:  Specimen refractive RI. Defaults to 1.47
            ti0: Working distance (immersion medium thickness) design value.
                 Defaults to 150
            tg:  Coverslip thickness experimental value. Defaults to 170.
            tg0: Coverslip thickness design value. Defaults to 170.
            sf (int):  oversampling factor to approximate pixel integration. Defaults to 3
            mode (int): if 0, returns oversampled PSF. Defaults to 1

    Returns:
        np.ndarray: The PSF volume.
"""

vectorialXYZFocalScan.__doc__ += _otherapidoc
scalarXYZFocalScan.__doc__ += _otherapidoc

__all__ = [
    "vectorial_psf",
    "scalar_psf",
    "vectorial_psf_centered",
    "scalar_psf_centered",
    "vectorialXYZFocalScan",
    "scalarXYZFocalScan",
    "normalize_params",
]
