import warnings
from typing import Sequence, Union

import _psfmodels
import numpy as np
from typing_extensions import Literal


def make_psf(
    z: Union[int, Sequence[float]] = 51,
    nx: int = 51,
    *,
    dxy: float = 0.05,
    dz: float = 0.05,
    pz: float = 0.0,
    NA: float = 1.4,
    wvl: float = 0.6,
    ns: float = 1.47,
    ni: float = 1.515,
    ni0: float = 1.515,
    tg: float = 170,
    tg0: float = 170,
    ng: float = 1.515,
    ng0: float = 1.515,
    ti0: float = 150.0,
    oversample_factor: int = 3,
    normalize: bool = True,
    model: Literal["vectorial", "scalar", "gaussian"] = "vectorial",
):
    """Compute microscope PSF.

    Select the PSF model using the `model` keyword argument. Can be one of:
        vectorial:  Vectorial PSF described by Aguet et al (2009).
        scalar:     Scalar PSF model described by Gibson and Lanni.
        gaussian:   Simple gaussian approximation.

    Parameters
    ----------
    z : Union[int, Sequence[float]]
        If an integer, z is interepreted as the number of z planes to calculate, and
        the point source always resides in the center of the z volume (at plane ~z//2).
        If a sequence (list, tuple, np.array), z is interpreted as a vector of Z
        positions at which the PSF is calculated (in microns, relative to
        coverslip).
        When an integer is provided, `dz` may be used to change the step size.
        If a sequence is provided, `dz` is ignored, since the sequence already implies
        absolute positions relative to the coverslip.
    nx : int
        XY size of output PSF in pixels, prefer odd numbers.
    dxy : float
        pixel size in sample space (microns)
    dz : float
        axial size in sample space (microns). Only used when `z` is an integer.
    pz : float
        point source z position above the coverslip, in microns.
    NA : float
        numerical aperture of the objective lens
    wvl : float
        emission wavelength (microns)
    ns : float
        sample refractive index
    ni : float
        immersion medium refractive index, experimental value
    ni0 : float
        immersion medium refractive index, design value
    tg : float
        coverslip thickness, experimental value (microns)
    tg0 : float
        coverslip thickness, design value (microns)
    ng : float
        coverslip refractive index, experimental value
    ng0 : float
        coverslip refractive index, design value
    ti0 : float
        working distance of the objective (microns)
    oversample_factor : int, optional
        oversampling factor to approximate pixel integration, by default 3
    normalize : bool
        Whether to normalize the max value to 1. By default, True.
    model : str
        PSF model to use.  Must be one of 'vectorial', 'scalar', 'gaussian'.
        By default 'vectorial'.

    Returns
    -------
    psf : np.ndarray
        The PSF array with dtype np.float64 and shape (len(zv), nx, nx)
    """
    if isinstance(z, (int, float)):
        zv: np.ndarray = _centered_zv(z, dz, pz)
    else:
        if dz != 0.05:
            warnings.warn("dz is ignored when providing a sequence for `z`.")
        zv = np.asarray(z)

    _args = set(_VALID_ARGS).difference({"zv", "nx", "dxy", "pz", "wvl"})
    kwargs = {k: v for k, v in locals().items() if k in _args}
    kwargs["sf"] = oversample_factor
    kwargs["mode"] = 1 if oversample_factor else 0

    if model == "vectorial":
        f = vectorial_psf
    elif model == "scalar":
        f = scalar_psf
    elif model == "gaussian":
        f = gaussian_psf
    else:
        raise ValueError(f"Unrecognized psf model: {model!r}")

    return f(zv=zv, nx=nx, dxy=dxy, pz=pz, wvl=wvl, params=kwargs, normalize=normalize)


# --------------------------------------------------------------------------------

_DEFAULT_PARAMS = {
    "NA": 1.4,  # numerical aperture
    "ng0": 1.515,  # coverslip RI design value
    "ng": 1.515,  # coverslip RI experimental value
    "ni0": 1.515,  # immersion medium RI design value
    "ni": 1.515,  # immersion medium RI experimental value
    "ns": 1.47,  # specimen refractive index (RI)
    "ti0": 150.0,  # microns, working distance (immersion medium thickness) design value
    "tg": 170.0,  # microns, coverslip thickness experimental value
    "tg0": 170.0,  # microns, coverslip thickness design value
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


def _normalize_params(mp):
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
            warnings.warn(
                f"parameter {key} is not one of the recognized keywords "
                "and is being ignored"
            )

    if _mp["NA"] >= _mp["ni"]:
        raise ValueError("NA must not be greater than immersion medium RI (ni).")
    return {key: val for key, val in _mp.items() if key in _VALID_ARGS}


def _validate_args(zv, dxy, pz):
    """Sanity checks for various arguments."""
    if isinstance(zv, (float, int)):
        zv = np.array([zv])
    elif isinstance(zv, (list, tuple)):
        zv = np.array(zv)
    elif not isinstance(zv, np.ndarray):
        raise ValueError("zd must be a scalar, iterable, or numpy array")
    if dxy <= 0:
        raise ValueError("dxy must be greater than 0")
    if pz < 0:
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
        nx (int, optional): XY size of output PSF in pixels, must be odd. Defaults to
        31. dxy (float, optional): XY Pixel size in sample space (microns). Defaults to
        0.05. pz (float, optional): Depth of point source relative to coverslip in
        microns.
                              Defaults to 0.
        wvl (float, optional): Emission wavelength in microns. Defaults to 0.6. params
        (dict, optional): Microscope parameters dict. See optional keys below. normalize
        (bool, optional): Normalize PSF peak to 1. Defaults to True.

        valid params (all floats unless stated, all distances in microns):
            NA:  Numerical Aperture. Defaults to 1.4 ng0: Coverslip RI design value.
            Defaults to 1.515 ng:  Coverslip RI experimental value. Defaults to 1.515
            ni0: Immersion medium RI design value. Defaults to 1.515 ni:  Immersion
            medium RI experimental value. Defaults to 1.515 ns:  Specimen refractive RI.
            Defaults to 1.47 ti0: Working distance (immersion medium thickness) design
            value.
                 Defaults to 150
            tg:  Coverslip thickness experimental value. Defaults to 170. tg0: Coverslip
            thickness design value. Defaults to 170. sf (int):  oversampling factor to
            approximate pixel integration. Defaults to 3 mode (int): if 0, returns
            oversampled PSF. Defaults to 1

    Returns:
        np.ndarray: The PSF volume.
 """


_DEFAULT_PARAMS = {
    "NA": 1.4,  # numerical aperture
    "ng0": 1.515,  # coverslip RI design value
    "ng": 1.515,  # coverslip RI experimental value
    "ni0": 1.515,  # immersion medium RI design value
    "ni": 1.515,  # immersion medium RI experimental value
    "ns": 1.47,  # specimen refractive index (RI)
    "ti0": 150.0,  # microns, working distance (immersion medium thickness) design value
    "tg": 170.0,  # microns, coverslip thickness experimental value
    "tg0": 170.0,  # microns, coverslip thickness design value
}


def vectorial_psf(zv=0, nx=31, dxy=0.05, pz=0.0, wvl=0.6, params=None, normalize=True):
    """Compute a vectorial model of the microscope point spread function."""
    zv = _validate_args(zv, dxy, pz)
    params = _normalize_params(params)
    _psf = _psfmodels.vectorial_psf(zv.copy(), int(nx), pz, wvl=wvl, dxy=dxy, **params)
    if normalize:
        _psf /= np.max(_psf)
    return _psf


def scalar_psf(zv=0, nx=31, dxy=0.05, pz=0, wvl=0.6, params=None, normalize=True):
    """Compute the scalar PSF model described by Gibson and Lanni."""
    zv = _validate_args(zv, dxy, pz)
    params = _normalize_params(params)
    _psf = _psfmodels.scalar_psf(zv.copy(), int(nx), pz, wvl=wvl, dxy=dxy, **params)
    if normalize:
        _psf /= np.max(_psf)
    return _psf


def gaussian_psf(zv=0, nx=31, dxy=0.05, pz=0, wvl=0.6, params=None, normalize=True):
    """Approximate 3D PSF as a gaussian.

    Parameters derived from Zhang et al (2007). https://doi.org/10.1364/AO.46.001819
    Using the paraxial approximation for NA < 0.7 and the Nonparaxial approximation
    for NA >= 0.7.
    """
    from scipy.stats import multivariate_normal

    if pz != 0:
        warnings.warn("pz != 0 currently does nothing for the gaussian approximation.")

    zv = _validate_args(zv, dxy, pz)
    params = _normalize_params(params)

    na = params["NA"]
    nimm = params["ni"]
    alpha = np.arcsin(na / nimm)
    cosa = np.cos(alpha)
    Kem = 2 * np.pi / wvl

    if na < 0.7:
        # paraxial
        sigma_xy = np.sqrt(2) / (Kem * na)
        sigma_z = (2 * np.sqrt(6) * nimm) / (Kem * na**2)

    else:
        # non-parax
        sigma_xy = (4 - 7 * cosa**1.5 + 3 * cosa**3.5) / (7 * (1 - cosa**1.5))
        sigma_xy = (1 / (nimm * Kem)) * sigma_xy**-0.5

        d = 4 * cosa**5 - 25 * cosa**3.5 + 42 * cosa**2.5 - 25 * cosa**1.5 + 4
        d = np.sqrt(6) * nimm * Kem * np.sqrt(d)
        sigma_z = 5 * np.sqrt(7) * (1 - cosa**1.5) / d

    _xycoords = _centered_zv(nx, dxy)
    y, z, x = np.meshgrid(_xycoords, zv, _xycoords)
    coords = np.column_stack([z.flat, y.flat, x.flat])

    sigma = np.array([sigma_z, sigma_xy, sigma_xy])
    _psf = multivariate_normal.pdf(coords, mean=[0, 0, 0], cov=np.diag(sigma**2))
    _psf = _psf.reshape((len(zv), nx, nx))

    if normalize:
        _psf /= np.max(_psf)
    return _psf


def vectorial_psf_deriv(
    zv=0, nx=31, dxy=0.05, pz=0.0, wvl=0.6, params=None, normalize=True
):
    """Compute a vectorial model of the microscope point spread function.

    also returns derivatives in dx, dy, dz.

    Returns:
        4-tuple of np.ndarrays: (_psf, dxp, dyp, dzp)

    """
    zv = _validate_args(zv, dxy, pz)
    params = _normalize_params(params)
    pixdxp = np.zeros((len(zv), nx, nx))
    pixdyp = np.zeros((len(zv), nx, nx))
    pixdzp = np.zeros((len(zv), nx, nx))
    _psf = _psfmodels.vectorial_psf_deriv(
        pixdxp, pixdyp, pixdzp, zv.copy(), int(nx), pz, wvl=wvl, dxy=dxy, **params
    )
    if normalize:
        _psf /= np.max(_psf)
    return _psf, pixdxp, pixdyp, pixdzp


def vectorial_psf_centered(nz, dz=0.05, **kwargs):
    """Compute a vectorial model of the microscope point spread function.

    The point source is always in the center of the output volume."""
    zv = _centered_zv(nz, dz, kwargs.get("pz", 0))
    return vectorial_psf(zv, **kwargs)


def scalar_psf_centered(nz, dz=0.05, **kwargs):
    """Compute the scalar PSF model described by Gibson and Lanni.

    The point source is always in the center of the output volume."""
    zv = _centered_zv(nz, dz, kwargs.get("pz", 0))
    return scalar_psf(zv, **kwargs)


def _centered_zv(nz, dz, pz=0) -> np.ndarray:
    lim = (nz - 1) * dz / 2
    return np.linspace(-lim + pz, lim + pz, nz)


vectorial_psf.__doc__ += _docstring + _paramdocs  # type: ignore
scalar_psf.__doc__ += _docstring + _paramdocs  # type: ignore
vectorial_psf_centered.__doc__ += _centerdocstring + _paramdocs  # type: ignore
scalar_psf_centered.__doc__ += _centerdocstring + _paramdocs  # type: ignore


def vectorialXYZFocalScan(mp, dxy, xy_size, zv, normalize=True, pz=0.0, wvl=0.6, zd=0):
    """Compute a vectorial model of the microscope point spread function.

    This function is merely here as a convenience to mimic the MicroscPSF-Py API.
    """
    return vectorial_psf(zv, xy_size, dxy, pz, wvl, mp, normalize)


def scalarXYZFocalScan(mp, dxy, xy_size, zv, normalize=True, pz=0.0, wvl=0.6, zd=0):
    """Compute the scalar PSF model described by Gibson and Lanni.

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
            sf (int):  oversampling factor to approximate pixel integration.
                Defaults to 3
            mode (int): if 0, returns oversampled PSF. Defaults to 1

    Returns:
        np.ndarray: The PSF volume.
"""

vectorialXYZFocalScan.__doc__ += _otherapidoc  # type: ignore
scalarXYZFocalScan.__doc__ += _otherapidoc  # type: ignore


def tot_psf(
    nx=127,
    nz=127,
    dxy=0.05,
    dz=0.05,
    pz=0,
    z_offset=0,
    x_offset=0,
    ex_wvl=0.488,
    em_wvl=0.525,
    ex_params=None,
    em_params=None,
    psf_func="vectorial",
):
    """Simlulate a total system psf with orthogonal illumination & detection.

    (e.g. SPIM)

    Args:
        nx (int, optional): XY size of output PSF in pixels, must be odd. Defaults to
        127. nz (int): Z size of output PSF in pixels, must be odd. Defaults to 127. dxy
        (float, optional): XY Pixel size in sample space (microns). Defaults to 0.05. dz
        (float, optional): Z step size of PSF in sample space. Defaults to 0.05 pz (int,
        optional): Depth of point source relative to coverslip in microns.
                            Defaults to 0.
        z_offset (int, optional): Defocus between the axial position of the excitation
            and the detection plane, with respect to the detection lens. Defaults to 0.
        x_offset (int, optional): Mismatch between the focal point of the excitation
        beam
            and the point source, along the propogation direction of the excitation
            beam. Defaults to 0.
        ex_wvl (float, optional): Emission wavelength in microns. Defaults to 0.488.
        em_wvl (float, optional): Emission wavelength in microns. Defaults to 0.525.
        ex_params ([type], optional): Excitation lens parameters dict. See keys below.
        em_params ([type], optional): Emission lens parameters dict. See keys below.
        psf_func (str, optional): The psf model to use.  Can be any of
            {'vectorial', 'scalar', or 'microscpsf'}.  Where 'microscpsf' uses the
            `gLXYZFocalScan` function from MicroscPSF-Py (if installed). Defaults to
            "vectorial".

        valid params (all floats unless stated, all distances in microns):
            NA:  Numerical Aperture. Defaults to 0.4 for excitation and 1.1 for emission
            ni0: Immersion medium RI design value. Defaults to 1.33 ni:  Immersion
            medium RI experimental value. Defaults to 1.33 ns:  Specimen refractive RI.
            Defaults to 1.33 tg:  Coverslip thickness experimental value. Defaults to 0
            (water immersion) tg0: Coverslip thickness design value. Defaults to 0
            (water immersion) ti0: Working distance (immersion medium thickness) design
            value.
                 Defaults to 150
            ng0: Coverslip RI design value. Defaults to 1.515 ng:  Coverslip RI
            experimental value. Defaults to 1.515

    Raises:
        ImportError: If `psf_func` == 'microscpsf' and MicroscPSF-Py cannot be imported
        ValueError: If `psf_func` is not one of {'vectorial', 'scalar', or 'microscpsf'}

    Returns:
        3-tuple of np.ndarrays:  ex_psf, em_psf, total_system_psf
    """

    _x_params = _DEFAULT_PARAMS.copy()
    _x_params.update({"ni0": 1.33, "ni": 1.33, "ns": 1.33, "tg": 0, "tg0": 0})
    _x_params["NA"] = 0.4
    _m_params = _x_params.copy()
    _m_params["NA"] = 1.1
    if ex_params is not None:
        _x_params.update(ex_params)
    if em_params is not None:
        _m_params.update(em_params)

    if psf_func.lower().startswith("microsc"):
        try:
            import microscPSF.microscPSF as msPSF
        except ImportError as e:
            raise ImportError(
                "Could not import MicroscPSF-Py. "
                'Install with "pip install MicroscPSF-Py"'
            ) from e

        f = msPSF.gLXYZFocalScan
        _x_params["zd0"] = _x_params.get("zd0", 200.0 * 1.0e3)
        _x_params["M"] = _x_params.get("M", 100)
        _m_params["zd0"] = _m_params.get("zd0", 200.0 * 1.0e3)
        _m_params["M"] = _m_params.get("M", 100)
    elif psf_func.lower().startswith("vectorial"):
        f = vectorialXYZFocalScan
    elif psf_func.lower().startswith("scalar"):
        f = scalarXYZFocalScan
    else:
        raise ValueError(
            'psf_func must be one of {"vectorial", "scalar", or "microscpsf"}'
        )

    lim = (nx - 1) * dxy / 2
    emzvec = np.linspace(-lim + x_offset, lim + x_offset, nx)
    if np.mod(z_offset / dz, 1) != 0:
        z_offset = dz * np.round(z_offset / dz)
        warnings.warn(
            "Not Implemented: z_offset must be an even multiple of dz. "
            "Coercing z_offset to nearest dz multiple: %s" % z_offset
        )
    _zoff = int(np.ceil(z_offset / dz))
    ex_nx = nz + 2 * np.abs(_zoff)
    exzvec = _centered_zv(nz, dz, pz)

    ex_psf = f(_x_params, dz, ex_nx, emzvec, pz=0, wvl=ex_wvl).T.sum(0)
    ex_psf = ex_psf[:nz] if _zoff >= 0 else ex_psf[-nz:]
    em_psf = f(_m_params, dxy, nx, exzvec, pz=pz, wvl=em_wvl)

    combined = ex_psf[:, :, np.newaxis] * em_psf
    return (ex_psf, em_psf, combined)


__all__ = [
    "vectorial_psf",
    "vectorial_psf_deriv",
    "scalar_psf",
    "vectorial_psf_centered",
    "scalar_psf_centered",
    "vectorialXYZFocalScan",
    "scalarXYZFocalScan",
    "tot_psf",
]
