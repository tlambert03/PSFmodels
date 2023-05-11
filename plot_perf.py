# benchmark the different implementations with perfplot
import contextlib

import numpy as np
import perfplot
from psfmodels.v2 import _aguet2009, _li2017

kernels = [
    lambda a: _li2017.radial_psf_li2017(a, a),
    lambda a: _aguet2009.radial_psf_aguet2009(a, a),
]
labels = ["numpy_li", "numpy_aguet"]

with contextlib.suppress(ImportError):
    from psfmodels.v2 import jax

    kernels.extend(
        (
            lambda a: jax.radial_psf_li2017(a, a),
            lambda a: jax.radial_psf_aguet2009(a, a),
        )
    )
    labels.extend(["jax_li", "jax_aguet"])

with contextlib.suppress(ImportError):
    from psfmodels.v2 import cupy

    kernels.extend(
        (
            lambda a: cupy.radial_psf_li2017(a, a),
            lambda a: cupy.radial_psf_aguet2009(a, a),
        )
    )
    labels.extend(["cupy_li", "cupy_aguet"])

perfplot.show(
    setup=lambda n: np.linspace(-3, 3, n),
    kernels=kernels,
    labels=labels,
    n_range=[2**k for k in range(6, 14)],
    xlabel="array size (n**2)",
    target_time_per_measurement=3,  # needs to be higher to allow jax to compile
    equality_check=None,
)
