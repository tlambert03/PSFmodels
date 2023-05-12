from __future__ import annotations

import inspect

import numpy
import torch
from torch.special import bessel_j0, bessel_j1


class ss:
    j0 = bessel_j0
    j1 = bessel_j1


from . import _aguet2009, _li2017

# recompile function with cupy as xp
_ns: dict = {}

_source = inspect.getsource(_li2017.radial_psf_li2017)
exec(compile(_source, __file__, "exec"), {"xp": torch, "ss": ss}, _ns)
radial_psf_li2017 = _ns["radial_psf_li2017"]

_source = inspect.getsource(_aguet2009._simpson)
exec(
    compile(_source, __file__, "exec"),
    {
        "xp": torch,
        "numpy": numpy,
        "ss": ss,
        "_simp_like": _aguet2009._simp_like,
    },
    _ns,
)
_simpson = _ns["_simpson"]

_source = inspect.getsource(_aguet2009.radial_psf_aguet2009)
exec(
    compile(_source, __file__, "exec"),
    {"xp": torch, "ss": ss, "_simpson": _simpson},
    _ns,
)
radial_psf_aguet2009 = _ns["radial_psf_aguet2009"]
