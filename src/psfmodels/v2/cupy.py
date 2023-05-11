from __future__ import annotations

import inspect

import cupy
import cupyx.scipy.special

from . import _aguet2009, _li2017

# recompile function with cupy as xp
_ns: dict = {}

_source = inspect.getsource(_li2017.radial_psf_li2017)
exec(compile(_source, __file__, "exec"), {"xp": cupy, "ss": cupyx.scipy.special}, _ns)
radial_psf_li2017 = _ns["radial_psf_li2017"]

_source = inspect.getsource(_aguet2009.radial_psf_aguet2009)
exec(compile(_source, __file__, "exec"), {"xp": cupy, "ss": cupyx.scipy.special}, _ns)
radial_psf_aguet2009 = _ns["radial_psf_aguet2009"]
