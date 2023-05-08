import numpy as np
import numpy.typing as npt


def nextpow2(n: int):
    """Return the next power of 2 greater than or equal to n."""
    return 2 ** int(np.ceil(np.log2(n)))


def bluestein_dft(
    x: npt.NDArray[np.complex_], f1: float, f2: float, fs: float, mout: int
):
    """Compute the Bluestein DFT of a signal.

    Parameters
    ----------
    x : np.ndarray
        Input signal. 2D complex array
    f1 : float
        Lower frequency bound
    f2 : float
        Upper frequency bound
    fs : float
        Sampling frequency
    mout : int
        Number of output samples
    """
    x = x.T
    n, m = x.shape
    f11 = f1 + (mout * fs + f2 - f1) / (2 * mout)
    f22 = f2 + (mout * fs + f2 - f1) / (2 * mout)
    a = np.exp(1j * 2 * np.pi * f11 / fs)
    w = np.exp(-1j * 2 * np.pi * (f22 - f11) / (mout * fs))
    h1 = np.arange(-m + 1, max(mout - 1, m - 1) + 1)
    mp = m + mout - 1
    h = w ** ((h1**2) / 2)
    ft = np.fft.fft(1 / h[:mp], nextpow2(mp))
    b1 = a ** (-np.arange(0, m)).reshape(-1, 1) * h[m - 1 : (2 * m) - 1].reshape(-1, 1)
    b = np.fft.fft(x * np.tile(b1.T, (n, 1)), nextpow2(mp))
    b = np.fft.ifft(b * np.tile(ft, (n, 1)))
    b = (b[:n, m - 1 : mp]) * np.tile(h[m - 1 : mp], (n, 1))
    l1 = np.linspace(0, mout - 1, mout)
    l = l1 / mout * (f22 - f11) + f11
    Mshift = np.tile(np.exp(-1j * 2 * np.pi * l * (-m / 2 + 1 / 2) / fs), (n, 1))
    return b * Mshift


def scalar_bluestein(
    gin: npt.NDArray[np.complex_],
    mxin: int,
    myin: int,
    pixelin: float,
    d: float,
    xstart: float,
    xend: float,
    ystart: float,
    yend: float,
    mxout: int,
    myout: int,
) -> npt.NDArray[np.complex_]:
    """Scalar diffraction computation method using Bluestein method

    Parameters
    ----------
    gin : npt.NDArray[np.complex_]
        Complex amplitude of the incident light beam
    mxin : int
        resolution of the input plane in transverse direction
    myin : int
        resolution of the input plane in longitudinal direction
    pixelin : float
        pixel size of the incident light beam
    d : float
        diffraction distance
    xstart : int
        computation range in transverse direction
    xend : int
        computation range in transverse direction
    ystart : int
        computation range in longitudinal direction
    yend : int
        computation range in longitudinal direction
    mxout : int
        resolution of the output plane in transverse direction
    myout : int
        resolution of the output plane in longitudinal direction

    Returns
    -------
    gout : npt.NDArray[np.complex_]
        Complex amplitude of the outgoing light beam
    """

    L0 = (myin - 1) * pixelin
    x0 = np.linspace(-L0 / 2, L0 / 2, mxin)
    y0 = np.linspace(-L0 / 2, L0 / 2, myin)
    [x0, y0] = np.meshgrid(x0, y0)

    # sometimes Nan in matlab code... appears unused.
    # pixelout = (xend - xstart) / (mxout - 1)

    x1 = np.linspace(xstart, xend, mxout)
    y1 = np.linspace(ystart, yend, myout)
    [x1, y1] = np.meshgrid(x1, y1)

    F0 = (
        np.exp(1j * k * d)
        / (1j * lamda * d)
        * np.exp(1j * k / 2 / d * (x1**2 + y1**2))
    )
    F = np.exp(1j * k / 2 / d * (x0**2 + y0**2))
    gout = gin * F

    fs = lamda * d / pixelin
    fy1 = ystart + fs / 2
    fy2 = yend + fs / 2
    gout = bluestein_dft(gout, fy1, fy2, fs, myout)

    fx1 = xstart + fs / 2
    fx2 = xend + fs / 2
    gout = bluestein_dft(gout, fx1, fx2, fs, mxout)
    gout = F0 * gout

    return gout


def scalar_cross_section():
    lamda = 550e-3
    k = 2 * np.pi / lamda

    # % Set parameters for aperture size and resolution
    my0 = 500
    mx0 = my0
    pixel0 = 30

    # % Create meshgrid of coordinates for aperture
    x = np.arange(-(mx0 - 1) / 2, (mx0 - 1) / 2 + 1)
    y = np.arange(-(my0 - 1) / 2, (my0 - 1) / 2 + 1)
    [xx, yy] = np.meshgrid(x, y)

    # % Scale coordinates to physical units (micrometers)
    x0 = xx * pixel0
    y0 = yy * pixel0

    # % Define aperture function as binary matrix using sign function
    aperture = np.sign(1 - np.sign(xx**2 + yy**2 - ((my0 - 1) / 2) ** 2))
    # % Set the initial wavefront to the aperture function
    A0 = aperture

    # % Propagate the wavefront to the observation plane using scalar diffraction
    f = 600e3
    g = A0 * np.exp(-1j * k / 2 / f * (x0**2 + y0**2))

    # % Define parameters for the observation plane
    L = 0.2e3
    ly = 201
    zstart = 571e3
    zfinish = 629e3

    # % Calculate the number of steps in propagation direction needed to reach the observation plane
    num = round((zfinish - zstart) / lamda + 1)
    num = round(num / ly)

    # % Adjust the start and finish points for the observation plane based on the number of steps
    zstart = 600e3 - (ly - 1) * num * lamda / 2
    zfinish = 600e3 + (ly - 1) * num * lamda / 2

    # % Define the x and y limits for the observation plane
    x1start = 0
    x1end = 0
    y1start = -L / 2
    y1end = L / 2

    # % Set the resolution of the observation plane
    mx1 = 1
    my1 = 1081

    # % Create meshgrid of coordinates for the observation plane
    x1 = np.linspace(x1start, x1end, mx1)
    y1 = np.linspace(y1start, y1end, my1)
    [x1, y1] = np.meshgrid(x1, y1)

    I_zs = []
    P_zs = []
    for d in np.linspace(zstart, zfinish, ly):
        # % Use Bluestein's algorithm to propagate the wavefront to the observation plane
        g1 = scalar_bluestein(
            g, mx0, my0, pixel0, d, x1start, x1end, y1start, y1end, mx1, my1
        )

        # % Calculate the intensity and phase of the wavefront at the observation plane
        I = np.abs(g1) ** 2
        P = np.angle(g1)

        # % Save the intensity and phase values along the y-axis for each propagation distance
        I_zs.append(I[:, mx1 // 2])
        P_zs.append(P[:, mx1 // 2])

    I_z = np.array(I_zs)
    P_z = np.array(P_zs)
    I_z = I_z / np.max(I_z)
    # [z1,y1]=np.meshgrid(
    #   zstart:((zfinish-zstart)/(ly-1)):zfinish,linspace(y1start,y1end,my1)
    # );

    return I_z, P_z
