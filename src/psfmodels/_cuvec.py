# import cupy as xp
# from cupyx.scipy.special import j0, j1
# from cupyx.scipy.ndimage import map_coordinates
from dataclasses import dataclass

import numpy as np
import numpy as xp
from scipy.ndimage import map_coordinates
from scipy.special import j0, j1


@dataclass
class Objective:
    na: float = 1.4  # numerical aperture
    coverslip_ri: float = 1.515  # coverslip RI experimental value (ng)
    coverslip_ri_spec: float = 1.515  # coverslip RI design value (ng0)
    immersion_medium_ri: float = 1.515  # immersion medium RI experimental value (ni)
    immersion_medium_ri_spec: float = 1.515  # immersion medium RI design value (ni0)
    specimen_ri: float = 1.47  # specimen refractive index (ns)
    working_distance: float = 150.0  # um, working distance, design value (ti0)
    coverslip_thickness: float = (
        170.0  # um, coverslip thickness experimental value (tg)
    )
    coverslip_thickness_spec: float = (
        170.0  # um, coverslip thickness design value (tg0)
    )

    @property
    def NA(self):
        return self.na

    @property
    def ng(self):
        return self.coverslip_ri

    @property
    def ng0(self):
        return self.coverslip_ri_spec

    @property
    def ni(self):
        return self.immersion_medium_ri

    @property
    def ni0(self):
        return self.immersion_medium_ri_spec

    @property
    def ns(self):
        return self.specimen_ri

    @property
    def ti0(self):
        return self.working_distance * 1e-6

    @property
    def tg(self):
        return self.coverslip_thickness * 1e-6

    @property
    def tg0(self):
        return self.coverslip_thickness_spec * 1e-6

    @property
    def half_angle(self):
        return np.arcsin(self.na / self.ni)


def wave_num(wavelength: float):
    # wavelength in microns
    2 * np.pi / (wavelength * 1e-6)


def _simp_like(arr):
    simp = xp.empty_like(arr)
    simp[::2] = 4
    simp[1::2] = 2
    simp[-1] = 1
    return simp


def simpson(
    p: Objective,
    theta: np.ndarray,
    constJ: np.ndarray,
    zv: np.ndarray,
    ci: float,
    zp: float,
    wave_num: float,
):
    sintheta = xp.sin(theta)
    costheta = xp.cos(theta)
    sqrtcostheta = xp.sqrt(costheta).astype("complex")
    ni2sin2theta = (p.ni**2 * sintheta**2).astype("complex")
    nsroot = xp.sqrt(p.ns**2 - ni2sin2theta)
    ngroot = xp.sqrt(p.ng**2 - ni2sin2theta)
    expW = xp.exp(
        1j
        * wave_num
        * (
            (ci - (zv[:, xp.newaxis, xp.newaxis] if zv.ndim else zv)) * p.ni * costheta
            + zp * nsroot
            + p.tg * ngroot
            - p.tg0 * xp.sqrt(p.ng0**2 - ni2sin2theta)
            - p.ti0 * xp.sqrt(p.ni0**2 - ni2sin2theta)
        )
    )

    simp = _simp_like(theta)

    ts1ts2 = (4.0 * p.ni * costheta * ngroot).astype("complex")
    tp1tp2 = ts1ts2.copy()
    tp1tp2 /= (p.ng * costheta + p.ni / p.ng * ngroot) * (
        p.ns / p.ng * ngroot + p.ng / p.ns * nsroot
    )
    ts1ts2 /= (p.ni * costheta + ngroot) * (ngroot + nsroot)

    # 2.0 factor: Simpson's rule
    bessel_0 = simp * j0(constJ[:, xp.newaxis] * sintheta) * sintheta * sqrtcostheta
    bessel_1 = simp * j1(constJ[:, xp.newaxis] * sintheta) * sintheta * sqrtcostheta

    with np.errstate(invalid="ignore"):
        bessel_2 = 2.0 * bessel_1 / (constJ[:, xp.newaxis] * sintheta) - bessel_0
    bessel_2[constJ == 0.0] = 0

    bessel_0 *= ts1ts2 + tp1tp2 / p.ns * nsroot
    bessel_1 *= tp1tp2 * p.ni / p.ns * sintheta
    bessel_2 *= ts1ts2 - tp1tp2 / p.ns * nsroot

    sum_I0 = xp.abs((expW * bessel_0).sum(-1))
    sum_I1 = xp.abs((expW * bessel_1).sum(-1))
    sum_I2 = xp.abs((expW * bessel_2).sum(-1))

    return xp.real(sum_I0**2 + 2.0 * sum_I1**2 + sum_I2**2)


def vectorial_rz(zv, nx=51, pos=(0, 0, 0), dxy=0.04, wvl=0.6, params={}, sf=3):
    p = Objective(**params)

    wave_num = 2 * np.pi / (wvl * 1e-6)

    xpos, ypos, zpos = pos

    # nz_ = len(z)
    xystep_ = dxy * 1e-6
    xymax = (nx * sf - 1) // 2

    # position in pixels
    xpos *= sf / xystep_
    ypos *= sf / xystep_
    rn = 1 + int(np.sqrt(xpos * xpos + ypos * ypos))
    rmax = int(np.ceil(np.sqrt(2.0) * xymax) + rn + 1)  # +1 for interpolation, dx, dy
    rvec = xp.arange(rmax) * xystep_ / sf
    constJ = wave_num * rvec * p.ni

    # CALCULATE
    # constant component of OPD
    ci = zpos * (1 - p.ni / p.ns) + p.ni * (p.tg0 / p.ng0 + p.ti0 / p.ni0 - p.tg / p.ng)

    nSamples = int(4 * (1.0 + p.half_angle * xp.max(constJ) / np.pi))
    nSamples = np.maximum(nSamples, 60)
    ud = 3.0 * sf

    theta = xp.arange(1, nSamples + 1) * p.half_angle / nSamples
    simpson_integral = simpson(p, theta, constJ, zv, ci, zpos, wave_num)
    step = p.half_angle / nSamples
    return 8.0 * np.pi / 3.0 * simpson_integral * (step / ud) ** 2

    # except xp.cuda.memory.OutOfMemoryError:
    #     integral = xp.empty((len(z), rmax))
    #     for k, zpos in enumerate(z):
    #         simp = simpson(p_, nSamples, constJ, zpos, ci, zp_)
    #         step = p.half_angle / nSamples
    #         integral[k] = 8.0 * np.pi / 3.0 * simp * (step / ud) ** 2
    #         del simp


def radius_map(shape, off=None):
    if off is not None:
        offy, offx = off
    else:
        off = (0, 0)
    ny, nx = shape
    yi, xi = xp.mgrid[:ny, :nx]
    yi = yi - (ny - 1) / 2 - offy
    xi = xi - (nx - 1) / 2 - offx
    return xp.hypot(yi, xi)


def rz_to_xyz(rz, xyshape, sf=3, off=None):
    """Use interpolation to create a 3D XYZ PSF from a 2D ZR PSF."""

    # Create XY grid of radius values.
    rmap = radius_map(xyshape, off) * sf
    nz = rz.shape[0]
    out = np.zeros((nz, *xyshape))
    for z in range(nz):
        o = map_coordinates(
            rz, xp.asarray([xp.ones(rmap.size) * z, rmap.ravel()]), order=1
        ).reshape(xyshape)
        out[z] = o.get() if hasattr(o, "get") else o
    return out


# def rz_to_xyz(rz, xyshape, sf=3, off=None):
#     """Use interpolation to create a 3D XYZ PSF from a 2D ZR PSF."""
#     # Create XY grid of radius values.
#     rmap = radius_map(xyshape, off) * sf
#     ny, nx = xyshape
#     nz, nr = rz.shape
#     ZZ, RR = xp.meshgrid(xp.arange(nz, dtype="float64"), rmap.ravel())
#     o = map_coordinates(rz, xp.array([ZZ.ravel(), RR.ravel()]))
#     return o.reshape((nx, ny, nz)).T


def vectorial_psf(
    zv, nx=51, ny=None, pos=(0, 0, 0), dxy=0.04, wvl=0.6, params={}, sf=3
):
    zv = xp.asarray(zv * 1e-6)  # convert to meters
    ny = ny or nx
    rz = vectorial_rz(zv, np.maximum(ny, nx), pos, dxy, wvl, params, sf)
    return rz_to_xyz(rz, (ny, nx), sf, off=np.array(pos[:2]) / (dxy * 1e-6))


if __name__ == "__main__":

    zv = np.linspace(-3, 3, 61)
    from time import perf_counter

    t0 = perf_counter()
    psf = vectorial_psf(zv, nx=512)
    t1 = perf_counter()
    print(psf.shape)
    print(t1 - t0)
    assert np.allclose(np.load("out.npy"), psf, atol=0.1)

    # import napari
    # napari.view_image(psf)
    # napari.run()
