"""Scalar and vectorial models of the microscope point spread function (PSF)."""

try:
    from ._version import version as __version__
except ImportError:
    __version__ = "unknown"
__author__ = "Talley Lambert"
__email__ = "talley.lambert@gmail.com"

from ._core import (
    confocal_psf,
    make_psf,
    scalar_psf,
    scalar_psf_centered,
    scalarXYZFocalScan,
    tot_psf,
    vectorial_psf,
    vectorial_psf_centered,
    vectorial_psf_deriv,
    vectorialXYZFocalScan,
)

__all__ = [
    "confocal_psf",
    "make_psf",
    "scalar_psf",
    "scalar_psf_centered",
    "scalarXYZFocalScan",
    "tot_psf",
    "vectorial_psf",
    "vectorial_psf_centered",
    "vectorial_psf_deriv",
    "vectorialXYZFocalScan",
]
