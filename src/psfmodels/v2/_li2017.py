from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as xp
import scipy.special as ss

if TYPE_CHECKING:
    import numpy.typing as npt


def radial_coords(dxy: float = 0.05, nxy: int = 64):
    """Return radial coordinates of the PSF array.

    Parameters
    ----------
    dxy : float
        Lateral pixel size, microns
    nxy : int
        Maximum lateral size of the PSF array, pixels
        This will be used to calculate the maximum radial coordinate.

    Returns
    -------
    np.ndarray
        Radial coordinates of the PSF array
    """
    rv_max = xp.sqrt(0.5 * nxy * nxy) + 1
    return xp.arange(0, rv_max * dxy, dxy)


def radial_psf_li2017(
    r: npt.ArrayLike | None = None,
    z: npt.ArrayLike | None = None,
    na: float = 1.4,
    wavelength: float = 0.6,
    *,
    pz: float = 0,
    ns: float = 1.47,
    ng0: float = 1.5,
    ng: float = 1.5,
    ni0: float = 1.5,
    ni: float = 1.5,
    ti0: float = 150,
    tg0: float = 170,
    tg: float = 170,
    min_wavelength: float = 0.436,
    num_basis: int = 100,
    rho_samples: int = 1000,
):
    """Compute a radial PSF slice using the Li, Xue, and Blu model.

    This implements a fast approximation  of the Gibson-Lanni  model, as described in:

    Li, J., Xue, F., & Blu, T. (2017). Fast and accurate three-dimensional point
    spread function computation for fluorescence microscopy. JOSA A, 34(6), 1029-1034

    Code modified from Kyle Douglass' implementation:

    https://github.com/kmdouglass/kmdouglass.github.io/blob/master/posts/implementing-a-fast-gibson-lanni-psf-solver-in-python/index.ipynb

    Parameters
    ----------
    r : np.ndarray
        Radial coordinates of the PSF array
    z : np.ndarray
        Axial coordinates of the PSF array
    NA : float, optional
        Numerical aperture, by default 1.4
    wavelength : float, optional
        Emission wavelength, microns, by default 0.610
    pz : float, optional
        Particle distance from coverslip, microns, by default 0
    ns : float, optional
        Specimen refractive index (RI), by default 1.33
    ng0 : float, optional
        Coverslip RI design value, by default 1.5
    ng : float, optional
        Coverslip RI experimental value, by default 1.5
    ni0 : float, optional
        Immersion medium RI design value, by default 1.5
    ni : float, optional
        Immersion medium RI experimental value, by default 1.5
    ti0 : float, optional
        Working distance (immersion medium thickness) design value, by default 150
    tg0 : float, optional
        Coverslip thickness design value, by default 170
    tg : float, optional
        Coverslip thickness experimental value, by default 170
    min_wavelength : float, optional
        Minimum emission wavelength, microns, by default 0.436
    num_basis : int, optional
        Number of rescaled Bessels that approximate the phase function, by default 100
    rho_samples : int, optional
        Number of pupil samples along radial direction, by default 1000

    Returns
    -------
    PSF_rz : np.ndarray
        3D PS
    """
    r = xp.asarray(r if r is not None else xp.linspace(0, 2, 41))
    z = xp.asarray(z if z is not None else xp.linspace(-3, 3, 61))

    # Scaling factors for the Fourier-Bessel series expansion
    scaling = na * (3 * xp.arange(1, num_basis + 1) - 2) * min_wavelength / wavelength

    # Wave number of emitted light.
    k = 2.0 * xp.pi / wavelength

    # Radial coordinates, pupil space
    max_rho = xp.min(xp.asarray([na, ns, ni, ni0, ng, ng0])) / na
    rho = xp.linspace(0, max_rho, rho_samples)

    # Define the wavefront aberration
    na2_rho2 = na**2 * rho**2
    # OPD in the sample.
    OPDs = pz * xp.sqrt(ns**2 - na2_rho2)
    # OPD in the immersion medium.
    ti = z.reshape(-1, 1) + ti0
    OPDi = ti * xp.sqrt(ni**2 - na2_rho2) - ti0 * xp.sqrt(ni0**2 - na2_rho2)
    # OPD in the coverslip.
    OPDg = tg * xp.sqrt(ng**2 - na2_rho2) - tg0 * xp.sqrt(ng0**2 - na2_rho2)
    # OPD in camera position.
    # a = na * zd0 / mag  # Aperture radius at the back focal plane.
    # OPDt = a * a * (zd0 - zd) * rho * rho / (2.0 * zd0 * zd)
    opd_tot = k * (OPDs + OPDi + OPDg)

    # Sample the phase
    # Shape is (number of z samples by number of rho samples)
    # phase = xp.cos(opd_tot) + 1j * xp.sin(opd_tot)
    phase = xp.exp(1j * opd_tot)

    # Define the basis of Bessel functions
    # Shape is (number of basis functions by number of rho samples)
    _q = scaling.reshape(-1, 1) * rho
    J = xp.asarray(ss.j0(_q), dtype=phase.dtype)  # cast only needed for torch

    # Compute the approximation to the sampled pupil phase by finding the least squares
    # solution to the complex coefficients of the Fourier-Bessel expansion.
    # Shape of C is (number of basis functions by number of z samples).
    # Note the matrix transposes to get the dimensions correct.
    # NOTE: rcond is for scipy >= 1.14.0
    C = xp.linalg.lstsq(J.T, phase.T, rcond=None)[0]

    b = k * na * r.reshape(-1, 1)

    # See equation 5 in Li, Xue, and Blu
    scaled_rho = scaling * max_rho
    _b_rho = b * max_rho
    R = (
        ss.j1(scaled_rho) * ss.j0(_b_rho) * scaled_rho
        - ss.j0(scaled_rho) * ss.j1(_b_rho) * _b_rho
    )
    R = R / (scaling * scaling - b * b)

    # The transpose places the axial direction along the first dimension of the array
    return (xp.abs(xp.asarray(R, dtype=C.dtype) @ C) ** 2).T
