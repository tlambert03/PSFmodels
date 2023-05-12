# benchmark the different implementations with perfplot
import contextlib

import numpy as np
import perfplot
from psfmodels.v2 import _aguet2009, _li2017

kernels = {
    "li_numpy": lambda a: _li2017.radial_psf_li2017(a, a),
    "aguet_numpy": lambda a: _aguet2009.radial_psf_aguet2009(a, a),
}

with contextlib.suppress(ImportError):
    from psfmodels.v2 import jax

    kernels["li_jax"] = lambda a: jax.radial_psf_li2017(a, a)
    kernels["aguet_jax"] = lambda a: jax.radial_psf_aguet2009(a, a)

with contextlib.suppress(ImportError):
    from psfmodels.v2 import cupy

    kernels["li_cupy"] = lambda a: cupy.radial_psf_li2017(a, a)
    kernels["aguet_cupy"] = lambda a: cupy.radial_psf_aguet2009(a, a)

with contextlib.suppress(ImportError):
    from psfmodels.v2 import torch

    kernels["li_torch"] = lambda a: torch.radial_psf_li2017(a, a)
    kernels["aguet_torch"] = lambda a: torch.radial_psf_aguet2009(a, a)


l, k = zip(*(x for x in kernels.items() if x[0].startswith("aguet_")))

perfplot.show(
    setup=lambda n: np.linspace(-3, 3, n),
    kernels=k,
    labels=l,
    n_range=[2**k for k in range(4, 6)],
    xlabel="array size (n**2)",
    target_time_per_measurement=4,  # needs to be higher to allow jax to compile
    equality_check=lambda x, y: np.allclose(x, y, atol=1e-3),  # torch is the most diff.
    title="Aguet2009"
)



l, k = zip(*(x for x in kernels.items() if x[0].startswith("li_")))

perfplot.show(
    setup=lambda n: np.linspace(-3, 3, n),
    kernels=k,
    labels=l,
    n_range=[2**k for k in range(4, 12)],
    xlabel="array size (n**2)",
    target_time_per_measurement=4,  # needs to be higher to allow jax to compile
    equality_check=lambda x, y: np.allclose(x, y, atol=1e-3),  # torch is the most diff.
    title="Li2017"
)

