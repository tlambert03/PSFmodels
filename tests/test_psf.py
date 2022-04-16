import numpy as np
import psfmodels as psfm
import pytest


@pytest.mark.parametrize("model", ["vectorial", "scalar", "gaussian"])
@pytest.mark.parametrize("z", [np.linspace(-2, 2, 15), 15])
def test_make_psf(model, z):
    p = psfm.make_psf(z, nx=31, model=model)
    assert p.shape == (15, 31, 31)


def test_vectorial_psf():
    zvec = np.linspace(-1, 1, 5)
    p = psfm.vectorial_psf(zvec, nx=31)
    assert p.shape == (5, 31, 31)


def test_vectorial_psf_deriv():
    zvec = np.linspace(-1, 1, 5)
    results = psfm.vectorial_psf_deriv(zvec, nx=31)
    assert len(results) == 4
    assert all(x.shape == (5, 31, 31) for x in results)


def test_vectorial_psf_centered():
    p = psfm.vectorial_psf_centered(nx=31, nz=5)
    assert p.shape == (5, 31, 31)


def test_scalar_psf():
    zvec = np.linspace(-1, 1, 5)
    p = psfm.scalar_psf(zvec, nx=31)
    assert p.shape == (5, 31, 31)


def test_scalar_psf_centered():
    p = psfm.scalar_psf_centered(nx=31, nz=5)
    assert p.shape == (5, 31, 31)
