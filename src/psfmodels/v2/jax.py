from __future__ import annotations

import inspect

import jax

from .. import _jax_bessel
from . import _li2017

# recompile function with jax.numpy as xp
_ns: dict = {}
_source = inspect.getsource(_li2017.radial_psf_li2017)
exec(compile(_source, __file__, "exec"), {"xp": jax.numpy, "ss": _jax_bessel}, _ns)
radial_psf_li2017 = jax.jit(_ns["radial_psf_li2017"])
