from __future__ import annotations

import inspect

import jax
import jax.numpy as jnp
import numpy

from .. import _jax_bessel
from . import _aguet2009, _li2017

# recompile function with jax.numpy as xp
_ns: dict = {}
_source = inspect.getsource(_li2017.radial_psf_li2017)
exec(compile(_source, __file__, "exec"), {"xp": jax.numpy, "ss": _jax_bessel}, _ns)
radial_psf_li2017 = jax.jit(_ns["radial_psf_li2017"])


def _simp_like(arr: jnp.ndarray):
    simp = jnp.empty_like(arr)

    simp = simp.at[::2].set(4)
    simp = simp.at[1::2].set(2)
    simp = simp.at[-1].set(1)
    return simp


_source = inspect.getsource(_aguet2009._simpson)
exec(
    compile(_source, __file__, "exec"),
    {
        "xp": jax.numpy,
        "numpy": numpy,
        "ss": _jax_bessel,
        "_simp_like": _simp_like,
    },
    _ns,
)
_simpson = jax.jit(_ns["_simpson"])

_source = inspect.getsource(_aguet2009.radial_psf_aguet2009)
exec(
    compile(_source, __file__, "exec"),
    {
        "xp": jax.numpy,
        "ss": _jax_bessel,
        "_simpson": _simpson,
    },
    _ns,
)
radial_psf_aguet2009 = jax.jit(_ns["radial_psf_aguet2009"])
