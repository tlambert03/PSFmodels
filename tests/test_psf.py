import unittest
import psfmodels as psfm
import numpy as np


class DeviceTest(unittest.TestCase):
    # This is a simple test that just tries to load the module

    def test_vectorial_psf(self):
        zvec = np.linspace(-1, 1, 5)
        p = psfm.vectorial_psf(zvec, nx=31)
        self.assertEqual(p.shape, (5, 31, 31))

    def test_vectorial_psf_deriv(self):
        zvec = np.linspace(-1, 1, 5)
        p, a, b, c = psfm.vectorial_psf_deriv(zvec, nx=31)
        self.assertEqual(p.shape, (5, 31, 31))
        self.assertEqual(a.shape, (5, 31, 31))
        self.assertEqual(b.shape, (5, 31, 31))
        self.assertEqual(c.shape, (5, 31, 31))

    def test_vectorial_psf_centered(self):
        p = psfm.vectorial_psf_centered(nx=31, nz=5)
        self.assertEqual(p.shape, (5, 31, 31))

    def test_scalar_psf(self):
        zvec = np.linspace(-1, 1, 5)
        p = psfm.scalar_psf(zvec, nx=31)
        self.assertEqual(p.shape, (5, 31, 31))

    def test_scalar_psf_centered(self):
        p = psfm.scalar_psf_centered(nx=31, nz=5)
        self.assertEqual(p.shape, (5, 31, 31))
