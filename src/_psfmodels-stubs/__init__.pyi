"""Scalar and Vectorial PSF Models."""
import numpy
import numpy.typing as npt

__all__ = ["scalar_psf", "vectorial_psf", "vectorial_psf_deriv"]

def scalar_psf(
    zv: npt.NDArray[numpy.float64],
    nx: int,
    pz: float,
    ti0: float,
    ni0: float,
    ni: float,
    tg0: float,
    tg: float,
    ng0: float,
    ng: float,
    ns: float,
    wvl: float,
    NA: float,
    dxy: float,
    sf: int = 3,
    mode: int = 1,
) -> npt.NDArray[numpy.float64]:
    """Compute scalar PSF model described by Gibson and Lanni.

    The model is described in F. Aguet et al., Opt. Express 17(8), pp. 6829-6848, 2009
    For more information and implementation details, see F. Aguet, Ph.D Thesis, Swiss
    Federal Institute of Technology, Lausanne (EPFL), 2009

    C++ code by Francois Aguet, 2009. Python bindings by Talley Lambert, 2019.

    Parameters
    ----------
    zv : np.ndarray
        Vector of Z positions at which PSF is calculated (in microns, relative to
        coverslip)
    nx : int
        XY size of output PSF in pixels, must be odd.
    pz : float
        point source z position above the coverslip in microns.
    ti0 : float
        working distance of the objective (microns)
    ni0 : float
        immersion medium refractive index, design value
    ni : float
        immersion medium refractive index, experimental value
    tg0 : float
        coverslip thickness, design value (microns)
    tg : float
        coverslip thickness, experimental value (microns)
    ng0 : float
        coverslip refractive index, design value
    ng : float
        coverslip refractive index, experimental value
    ns : float
        sample refractive index
    wvl : float
        emission wavelength (microns)
    NA : float
        numerical aperture
    dxy : float
        pixel size in sample space (microns)
    sf : int, optional
        oversampling factor to approximate pixel integration, by default 3
    mode : int, optional
        if 0, returns oversampled PSF, by default 1

    Returns
    -------
    npt.NDArray[numpy.float64]
        PSF with type np.float64 and shape (len(zv), nx, nx)
    """

def vectorial_psf(
    zv: npt.NDArray[numpy.float64],
    nx: int,
    pz: float,
    ti0: float,
    ni0: float,
    ni: float,
    tg0: float,
    tg: float,
    ng0: float,
    ng: float,
    ns: float,
    wvl: float,
    NA: float,
    dxy: float,
    sf: int = 3,
    mode: int = 1,
) -> npt.NDArray[numpy.float64]:
    """Compute vectorial microscope point spread function model.

    The model is described in F. Aguet et al., Opt. Express 17(8), pp. 6829-6848, 2009
    For more information and implementation details, see F. Aguet, Ph.D Thesis, Swiss
    Federal Institute of Technology, Lausanne (EPFL), 2009

    C++ code by Francois Aguet, 2009. Python bindings by Talley Lambert, 2019.

    Parameters
    ----------
    zv : np.ndarray
        Vector of Z positions at which PSF is calculated (in microns, relative to
        coverslip)
    nx : int
        XY size of output PSF in pixels, must be odd.
    pz : float
        point source z position above the coverslip in microns.
    ti0 : float
        working distance of the objective (microns)
    ni0 : float
        immersion medium refractive index, design value
    ni : float
        immersion medium refractive index, experimental value
    tg0 : float
        coverslip thickness, design value (microns)
    tg : float
        coverslip thickness, experimental value (microns)
    ng0 : float
        coverslip refractive index, design value
    ng : float
        coverslip refractive index, experimental value
    ns : float
        sample refractive index
    wvl : float
        emission wavelength (microns)
    NA : float
        numerical aperture
    dxy : float
        pixel size in sample space (microns)
    sf : int, optional
        oversampling factor to approximate pixel integration, by default 3
    mode : int, optional
        if 0, returns oversampled PSF, by default 1

    Returns
    -------
    npt.NDArray[numpy.float64]
        PSF with type np.float64 and shape (len(zv), nx, nx)
    """

def vectorial_psf_deriv(
    pixdxp: npt.NDArray[numpy.float64],
    pixdyp: npt.NDArray[numpy.float64],
    pixdzp: npt.NDArray[numpy.float64],
    zv: npt.NDArray[numpy.float64],
    nx: int,
    pz: float,
    ti0: float,
    ni0: float,
    ni: float,
    tg0: float,
    tg: float,
    ng0: float,
    ng: float,
    ns: float,
    wvl: float,
    NA: float,
    dxy: float,
    sf: int = 3,
    mode: int = 1,
) -> npt.NDArray[numpy.float64]:
    """Compute vectorial point spread function model, and return derivatives.

    Parameters
    ----------
    pixdxp : npt.NDArray[numpy.float64]
        Derivative of pixel x position with respect to z position
    pixdyp : npt.NDArray[numpy.float64]
        Derivative of pixel y position with respect to z position
    pixdzp : npt.NDArray[numpy.float64]
        Derivative of pixel z position with respect to z position
    zv : np.ndarray
        Vector of Z positions at which PSF is calculated (in microns, relative to
        coverslip)
    nx : int
        XY size of output PSF in pixels, must be odd.
    pz : float
        point source z position above the coverslip in microns.
    ti0 : float
        working distance of the objective (microns)
    ni0 : float
        immersion medium refractive index, design value
    ni : float
        immersion medium refractive index, experimental value
    tg0 : float
        coverslip thickness, design value (microns)
    tg : float
        coverslip thickness, experimental value (microns)
    ng0 : float
        coverslip refractive index, design value
    ng : float
        coverslip refractive index, experimental value
    ns : float
        sample refractive index
    wvl : float
        emission wavelength (microns)
    NA : float
        numerical aperture
    dxy : float
        pixel size in sample space (microns)
    sf : int, optional
        oversampling factor to approximate pixel integration, by default 3
    mode : int, optional
        if 0, returns oversampled PSF, by default 1

    Returns
    -------
    numpy.ndarray
        PSF with type np.float64 and shape (len(zv), nx, nx)
    """
