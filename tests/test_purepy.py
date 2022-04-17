import numpy as np
import psfmodels as psfm
from psfmodels import _cuvec as pure


def test_equality():
    N = 101
    dx = 0.001
    zv = np.linspace(-2, 2, N)
    a = psfm.vectorial_psf(zv, nx=N, dxy=dx)
    b = pure.vectorial_psf(zv, nx=N, dxy=dx)
    # note, similarity gets worse as dx goes up...
    # hints at an interpolation difference
    np.testing.assert_allclose(a, b, rtol=0.005)
