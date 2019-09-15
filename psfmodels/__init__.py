try:
    import _psfmodels as library
except ImportError:
    raise ImportError("could not import precompiled library for psfmodels")
import numpy as np
import warnings


DEFAULT_PARAMS = {
    "NA": 1.4,  # numerical aperture
    "ng0": 1.515,  # coverslip RI design value
    "ng": 1.515,  # coverslip RI experimental value
    "ni0": 1.515,  # immersion medium RI design value
    "ni": 1.515,  # immersion medium RI experimental value
    "ns": 1.33,  # specimen refractive index (RI)
    "ti0": 150,  # microns, working distance (immersion medium thickness) design value
    "tg": 170,  # microns, coverslip thickness experimental value
    "tg0": 170,  # microns, coverslip thickness design value
}

VALID_ARGS = [
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
    _mp = DEFAULT_PARAMS.copy()
    if mp is not None:
        if isinstance(mp, dict):
            _mp.update(mp)
        else:
            raise ValueError("mp argument must be dict of microscope params")
    out = {}
    for key, val in _mp.items():
        if key in VALID_ARGS:
            out[key] = val
        else:
            warnings.warn(
                (
                    "parameter {} is not one of the recognized "
                    "keywords and is being ignored"
                ).format(key)
            )

    return {key: val for key, val in _mp.items() if key in VALID_ARGS}


def _validate_args(zv, dxy, pz):
    if isinstance(zv, (float, int)):
        zv = np.array([zv])
    elif isinstance(zv, (list, tuple)):
        zv = np.array(zv)
    elif not isinstance(zv, np.ndarray):
        raise ValueError("zd must be a scalar, iterable, or numpy array")
    if not dxy > 0:
        raise ValueError("dxy must be greater than 0")
    if not pz >= 0:
        raise ValueError("pz should be >= 0")
    return zv


def vectorial_psf(zv=0, nx=31, dxy=0.05, pz=0, wvl=0.6, mp=None, normalize=True):
    zv = _validate_args(zv, dxy, pz)
    mp = normalize_params(mp)
    _psf = library.vectorial_psf(zv, int(nx), pz, wvl=wvl, dxy=dxy, **mp)
    if normalize:
        _psf /= np.max(_psf)
    return _psf


def scalar_psf(zv=0, nx=31, dxy=0.05, pz=0, wvl=0.6, mp=None, normalize=True):
    zv = _validate_args(zv, dxy, pz)
    mp = normalize_params(mp)
    _psf = library.scalar_psf(zv, int(nx), pz, wvl=wvl, dxy=dxy, **mp)
    if normalize:
        _psf /= np.max(_psf)
    return _psf


def _centered_zv(nz, dz, pz):
    lim = (nz - 1) * dz / 2
    return np.linspace(-lim + pz, lim + pz, nz)


def vectorial_psf_centered(nz, dz=0.05, **kwargs):
    zv = _centered_zv(nz, dz, kwargs.get("pz", 0))
    return vectorial_psf(zv, **kwargs)


def scalar_psf_centered(nz, dz=0.05, **kwargs):
    zv = _centered_zv(nz, dz, kwargs.get("pz", 0))
    return scalar_psf(zv, **kwargs)


def vectorialXYZFocalScan(mp, dxy, xy_size, zv, normalize=True, pz=0.0, wvl=0.6, zd=0):
    return vectorial_psf(zv, xy_size, dxy, pz, wvl, mp, normalize)


def scalarXYZFocalScan(mp, dxy, xy_size, zv, normalize=True, pz=0.0, wvl=0.6, zd=0):
    return scalar_psf(zv, xy_size, dxy, pz, wvl, mp, normalize)
