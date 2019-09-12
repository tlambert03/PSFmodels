import unittest
from vectorialpsf import make_psf
import numpy as np


class DeviceTest(unittest.TestCase):
    # This is a simple test that just tries to load the module

    def test_make_psf(self):
        zvec = np.linspace(-1, 1, 5)
        p = make_psf(zvec, nx=31)
        print(p.shape)
