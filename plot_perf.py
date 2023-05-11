# benchmark the different implementations with perfplot
import numpy as np
import perfplot
from psfmodels.v2 import _li2017, jax

perfplot.show(
    setup=lambda n: np.linspace(-3, 3, n),
    kernels=[
        lambda a: _li2017.radial_psf_li2017(a, a),
        lambda a: jax.radial_psf_li2017(a, a),
    ],
    labels=["numpy", "jax"],
    n_range=[2**k for k in range(4, 13)],
    target_time_per_measurement=3,  # needs to be higher to allow jax to compile
)
