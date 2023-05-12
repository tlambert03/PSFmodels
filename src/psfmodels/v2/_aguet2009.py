import numpy
import numpy as xp
import scipy.special as ss


def _simp_like(arr):
    simp = xp.empty_like(arr)
    simp[::2] = 4
    simp[1::2] = 2
    simp[-1] = 1
    return simp


def radial_psf_aguet2009(
    r: xp.ndarray | None = None,
    z: xp.ndarray | None = None,
    na: float = 1.4,
    wavelength: float = 0.6,
    *,
    pz: float = 0,
    ns: float = 1.47,
    ng0: float = 1.5,
    ng: float = 1.5,
    ni0: float = 1.5,
    ni: float = 1.5,
    ti0: float = 150,
    tg0: float = 170,
    tg: float = 170,
):
    rvec = xp.asarray(r if r is not None else xp.linspace(0, 2, 41))
    z = xp.asarray(z if z is not None else xp.linspace(-3, 3, 61))

    rvec = rvec * 1e-6
    z = z * 1e-6
    wavelength *= 1e-6
    ti0 *= 1e-6
    tg0 *= 1e-6
    tg *= 1e-6

    # Wave number of emitted light.
    k = 2.0 * xp.pi / wavelength

    constJ: xp.ndarray = k * rvec * ni

    # constant component of OPD
    ci = pz * (1 - ni / ns) + ni * (tg0 / ng0 + ti0 / ni0 - tg / ng)

    half_angle = xp.arcsin(xp.asarray([na / ni]))  # array cast is for torch
    nSamples = xp.asarray(4 * (1.0 + half_angle * xp.max(constJ) / xp.pi), dtype=int)
    nSamples = xp.maximum(nSamples, xp.asarray([60]))
    nSamples = 60

    step = half_angle / nSamples
    theta = xp.arange(1, nSamples + 1) * step
    simpson_integral = _simpson(
        theta, constJ, z, ci, pz, k, ni, ns, ng, tg, tg0, ni0, ti0, ng0
    )
    return 8.0 * xp.pi / 3.0 * simpson_integral * step**2


def _simpson(
    theta: xp.ndarray,
    constJ: xp.ndarray,
    zv: xp.ndarray,
    ci: float,
    zp: float,
    wave_num: float,
    ni: float,
    ns: float,
    ng: float,
    tg: float,
    tg0: float,
    ni0: float,
    ti0: float,
    ng0: float,
):
    # L_theta calculation
    sintheta = xp.sin(theta)
    costheta = xp.cos(theta)
    sqrtcostheta = xp.asarray(xp.sqrt(costheta), dtype=xp.complex128)
    ni2sin2theta = ni**2 * sintheta**2
    nsroot = xp.sqrt(ns**2 - ni2sin2theta)
    ngroot = xp.sqrt(ng**2 - ni2sin2theta)
    _z = zv[:, None, None] if zv.ndim else zv
    L0 = (
        ni * (ci - _z) * costheta
        + zp * nsroot
        + tg * ngroot
        - tg0 * xp.sqrt(ng0**2 - ni2sin2theta)
        - ti0 * xp.sqrt(ni0**2 - ni2sin2theta)
    )
    expW = xp.exp(1j * wave_num * L0)

    simp = xp.asarray(_simp_like(theta))

    ts1ts2 = xp.asarray(4.0 * ni * costheta * ngroot, dtype=xp.complex128)
    tp1tp2 = ts1ts2.clone() if hasattr(ts1ts2, "clone") else ts1ts2.copy()  # for torch
    tp1tp2 /= (ng * costheta + ni / ng * ngroot) * (ns / ng * ngroot + ng / ns * nsroot)
    ts1ts2 /= (ni * costheta + ngroot) * (ngroot + nsroot)

    # 2.0 factor: Simpson's rule
    bessel_0 = simp * ss.j0(constJ[:, None] * sintheta) * sintheta * sqrtcostheta
    bessel_1 = simp * ss.j1(constJ[:, None] * sintheta) * sintheta * sqrtcostheta

    with numpy.errstate(invalid="ignore"):
        bessel_2 = 2.0 * bessel_1 / (constJ[:, None] * sintheta) - bessel_0

    bessel_2 = xp.where((constJ == 0.0)[:, None], 0, bessel_2)

    bessel_0 *= ts1ts2 + tp1tp2 / ns * nsroot
    bessel_1 *= tp1tp2 * ni / ns * sintheta
    bessel_2 *= ts1ts2 - tp1tp2 / ns * nsroot

    sum_I0 = xp.abs((expW * bessel_0).sum(-1))
    sum_I1 = xp.abs((expW * bessel_1).sum(-1))
    sum_I2 = xp.abs((expW * bessel_2).sum(-1))

    return xp.real(sum_I0**2 + 2.0 * sum_I1**2 + sum_I2**2)
