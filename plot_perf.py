# benchmark the different implementations with perfplot
import contextlib

import numpy as np
import perfplot
from psfmodels.v2 import _li2017

kernels = [lambda a: _li2017.radial_psf_li2017(a, a)]
labels = ["numpy"]

with contextlib.suppress(ImportError):
    from psfmodels.v2 import jax

    kernels.append(lambda a: jax.radial_psf_li2017(a, a))
    labels.append("jax")

with contextlib.suppress(ImportError):
    from psfmodels.v2 import cupy

    kernels.append(lambda a: cupy.radial_psf_li2017(a, a))
    labels.append("cupy")

perfplot.show(
    setup=lambda n: np.linspace(-3, 3, n),
    kernels=kernels,
    labels=labels,
    n_range=[2**k for k in range(6, 14)],
    xlabel="array size (n**2)",
    target_time_per_measurement=4,  # needs to be higher to allow jax to compile
)
