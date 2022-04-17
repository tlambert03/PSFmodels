from enum import Enum
from typing import TYPE_CHECKING

from typing_extensions import Annotated

from . import _core

if TYPE_CHECKING:
    import napari.types


class PSFModel(Enum):
    vectorial = "vectorial"
    scalar = "scalar"
    gaussian = "gaussian"


dRI = {"min": 1, "max": 1.7, "step": 0.001}


def make_psf(
    model: Annotated[PSFModel, {"label": "PSF Model"}] = PSFModel.vectorial,
    nz: Annotated[int, {"min": 1, "max": 1025, "label": "Num. pixels (Z)"}] = 101,
    nx: Annotated[int, {"min": 1, "max": 1025, "label": "Num. pixels (XY)"}] = 101,
    dz: Annotated[
        float, {"min": 0.01, "max": 0.5, "step": 0.01, "label": "Pix size Z (µm)"}
    ] = 0.05,
    dxy: Annotated[
        float, {"min": 0.01, "max": 0.5, "step": 0.01, "label": "Pix size XY (µm)"}
    ] = 0.05,
    NA: Annotated[float, {"min": 0.1, "max": 1.7, "step": 0.01}] = 1.4,
    wvl: Annotated[
        float, {"min": 0.3, "max": 0.9, "step": 0.005, "label": "Wavelength (µm)"}
    ] = 0.6,
    pz: Annotated[float, {"max": 200, "label": "Depth from CS (µm)"}] = 0.0,
    ns: Annotated[float, {**dRI, "label": "Sample RI"}] = 1.47,
    ni: Annotated[float, {**dRI, "label": "Imm. RI"}] = 1.515,
    ni0: Annotated[float, {**dRI, "label": "Imm. RI (spec)"}] = 1.515,
    tg: Annotated[int, {"min": 0, "max": 200, "label": "CS thickness (µm)"}] = 170,
    tg0: Annotated[int, {"min": 0, "max": 200, "label": "CS thickness (spec)"}] = 170,
    ng: Annotated[float, {**dRI, "label": "CS RI"}] = 1.515,
    ng0: Annotated[float, {**dRI, "label": "CS RI (spec)"}] = 1.515,
) -> "napari.types.ImageData":
    """Generate 3D microscope PSF

    Select from one of the following PSF models:
        vectorial:  Vectorial PSF (Aguet et al, 2009)
        scalar:     Scalar PSF model (Gibson & Lanni, 1992)
        gaussian:   Gaussian approximation (Zhang et al, 2007)

    Parameters
    ----------
    model : str
        PSF model to use
            vectorial - Vectorial PSF (Aguet et al, 2009)
            scalar - Scalar PSF (Gibson & Lanni, 1992)
            gaussian - Gaussian approximation (Zhang et al, 2007)
    nz : Union[int, Sequence[float]]
        Number of z planes
    nx : int
        Number of XY pixels, (prefer odd numbers)
    dxy : float
        pixel size (µm)
    dz : float
        Z step size (µm)
    pz : float
        Point source position above the coverslip (µm)
    NA : float
        Numerical aperture of the objective lens
    wvl : float
        Emission wavelength (µm)
    ns : float
        Sample refractive index
    ni : float
        Immersion medium refractive index
    ni0 : float
        Immersion medium refractive index (design value)
    tg : float
        Coverslip thickness (µm)
    tg0 : float
        Coverslip thickness (design value, µm)
    ng : float
        Coverslip refractive index
    ng0 : float
        Coverslip refractive index (design value)
    """
    kwargs = locals().copy()
    kwargs["model"] = kwargs["model"].value
    kwargs["z"] = kwargs.pop("nz")
    return _core.make_psf(**kwargs)
