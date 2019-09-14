#include "vectorialPSF.cpp"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

double defaultxp[3] = {0.0, 0.0, 0.0};

py::array_t<double> make_psf(py::array_t<double> zvec, int nx = 101,
                             double zdepth = 0, double ti0 = 190,
                             double ni0 = 1.518, double ni = 1.518,
                             double tg0 = 170, double tg = 170,
                             double ng0 = 1.515, double ng = 1.515,
                             double ns = 1.47, double wvl = 0.550,
                             double mag = 100, double na = 1.45,
                             double pixel = 6.5, int sf = 3, int mode = 1) {
                               
  py::buffer_info buf1 = zvec.request();
  double *ptr1 = (double *)buf1.ptr;
  if (buf1.ndim != 1)
    throw std::runtime_error("zvec must be a 1-dimensional array");

  // convert microns to meters
  for (size_t i = 0; i < zvec.size(); i++) {
    ptr1[i] *= 1e-6;
  }

  py::print(zvec);
  py::print(zdepth);
  double xp[] = {0.0, 0.0, zdepth * 1e-6};

  parameters p;

  p.ti0 = ti0 * 1e-6;
  p.ni0 = ni0;
  p.ni = ni;
  p.tg0 = tg0 * 1e-6;
  p.tg = tg * 1e-6;
  p.ng0 = ng0;
  p.ng = ng;
  p.ns = ns;
  p.lambda = wvl * 1e-6;
  p.M = mag;
  p.NA = na;
  p.pixelSize = pixel * 1e-6;
  p.sf = sf;
  p.mode = mode;
  p.k0 = 2 * PI / p.lambda;
  p.ni0_2 = p.ni0 * p.ni0;
  p.ni_2 = p.ni * p.ni;
  p.ng0_2 = p.ng0 * p.ng0;
  p.ng_2 = p.ng * p.ng;
  p.ns_2 = p.ns * p.ns;
  p.alpha = asin(p.NA / p.ni);
  p.NA_2 = p.NA * p.NA;

  int nz = zvec.shape(0);
  VectorialPSF psf = VectorialPSF(xp, ptr1, nz, nx, p);
  psf.calculatePSF();
  return py::array_t<double, py::array::c_style | py::array::forcecast>(
      std::vector<ptrdiff_t>{nz, nx, nx}, &psf.pixels_[0]);
}

py::array_t<double>
make_centered_psf(int nx = 101, int nz = 31, double range = 6,
                  double zdepth = 0, double ti0 = 190, double ni0 = 1.518,
                  double ni = 1.518, double tg0 = 170, double tg = 170,
                  double ng0 = 1.515, double ng = 1.515, double ns = 1.47,
                  double wvl = 0.550, double mag = 100, double na = 1.45,
                  double pixel = 6.5, int sf = 3, int mode = 1) {

  double _zdepth = zdepth;
  double lim = range / 2;
  std::vector<double> zvecv = linspace(-lim + _zdepth, lim + _zdepth, nz);
  py::array_t<double> zvec(std::vector<ptrdiff_t>{nz}, &zvecv[0]);

  return make_psf(zvec, nx, zdepth, ti0, ni0, ni, tg0, tg, ng0, ng, ns, wvl,
                  mag, na, pixel, sf, mode);
}

PYBIND11_MODULE(vectorialpsf, m) {
  m.doc() = "Vectorial PSF generation"; // optional module docstring

  m.def("make_psf", &make_psf, R"pbdoc(
    Computes a vectorial microscope point spread function model.
    
    The model is described in [1].
    For more information and implementation details, see [2].

    Args:
        zvec (np.ndarray): Vector of Z positions at which PSF is calculated (in microns, relative to coverslip) 
        nx (int): XY size of output PSF in pixels, must be odd.
        zdepth (float): depth of point source relative to coverslip (in microns).
        ti0 (float): (optional, default: 1.515)  working distance of the objective (microns)
        ni0 (float): (optional, default: 1.518) immersion medium refractive index, design value
        ni (float): (optional, default: 1.518) immersion medium refractive index, experimental value
        tg0 (float): (optional, default: 170) coverslip thickness, design value (microns)
        tg (float): (optional, default: 170) coverslip thickness, experimental value (microns)
        ng0 (float): (optional, default: 1.515) coverslip refractive index, design value
        ng (float): (optional, default: 1.515) coverslip refractive index, experimental value
        ns (float): (optional, default: 1.47) sample refractive index
        wvl (float): (optional, default: 0.55) emission wavelength (microns)
        mag (float): (optional, default: 100) magnification
        na (float): (optional, default: 1.45) numerical aperture
        pixel (float): (optional, default: 6.5) physical size (width) of the camera pixels (microns)
        sf (int): (optional, default: 3) oversampling factor to approximate pixel integration
        mode (int): (optional, default: 1) if 0, returns oversampled PSF

    Returns:
        np.ndarray: PSF with type np.float64 and shape (len(zvec), nx, nx)

    [1] F. Aguet et al., Opt. Express 17(8), pp. 6829-6848, 2009
    [2] F. Aguet, Ph.D Thesis, Swiss Federal Institute of Technology, Lausanne (EPFL), 2009
    C code by Francois Aguet, 2009. Python binding by Talley Lambert, 2019.
    )pbdoc",
        py::arg("zvec"), py::arg("nx") = 101, py::arg("zdepth") = 0,
        py::arg("ti0") = 190, py::arg("ni0") = 1.518, py::arg("ni") = 1.518,
        py::arg("tg0") = 170, py::arg("tg") = 170, py::arg("ng0") = 1.515,
        py::arg("ng") = 1.515, py::arg("ns") = 1.47, py::arg("wvl") = 0.550,
        py::arg("mag") = 100, py::arg("na") = 1.45, py::arg("pixel") = 6.5,
        py::arg("sf") = 3, py::arg("mode") = 1);

  m.def("make_centered_psf", &make_centered_psf, R"pbdoc(
    Computes a vectorial microscope point spread function model.
    
    this is a convenience interface for make_psf, where the point source is always
    assumed to be in the middle of a volume with size nz, and range `range`, 
    with the point source at `zdepth` microns from the coverslip.

    The model is described in [1].
    For more information and implementation details, see [2].

    Args:
        nx (int): XY size of output PSF in pixels, must be odd.
        nz (int): XZ size of output PSF in pixels, must be odd.
        range (float): physical size of Z range (in microns)
        zdepth (float): depth of point source relative to coverslip (in microns)
        ti0 (float): (optional, default: 1.515)  working distance of the objective (microns)
        ni0 (float): (optional, default: 1.518) immersion medium refractive index, design value
        ni (float): (optional, default: 1.518) immersion medium refractive index, experimental value
        tg0 (float): (optional, default: 170) coverslip thickness, design value (microns)
        tg (float): (optional, default: 170) coverslip thickness, experimental value (microns)
        ng0 (float): (optional, default: 1.515) coverslip refractive index, design value
        ng (float): (optional, default: 1.515) coverslip refractive index, experimental value
        ns (float): (optional, default: 1.47) sample refractive index
        wvl (float): (optional, default: 0.55) emission wavelength (microns)
        mag (float): (optional, default: 100) magnification
        na (float): (optional, default: 1.45) numerical aperture
        pixel (float): (optional, default: 6.5) physical size (width) of the camera pixels (microns)
        sf (int): (optional, default: 3) oversampling factor to approximate pixel integration
        mode (int): (optional, default: 1) if 0, returns oversampled PSF

    Returns:
        np.ndarray: PSF with type np.float64 and shape (nz, nx, nx)
                    particle will always be at the center of the output volume
                    (i.e. [(nz+1)//2, (nx+1)//2, (nx+1)//2])

    [1] F. Aguet et al., Opt. Express 17(8), pp. 6829-6848, 2009
    [2] F. Aguet, Ph.D Thesis, Swiss Federal Institute of Technology, Lausanne (EPFL), 2009
    C code by Francois Aguet, 2009. Python binding by Talley Lambert, 2019.
    )pbdoc",
        py::arg("nx") = 101, py::arg("nz") = 31, py::arg("range") = 6,
        py::arg("zdepth") = 0, py::arg("ti0") = 190, py::arg("ni0") = 1.518,
        py::arg("ni") = 1.518, py::arg("tg0") = 170, py::arg("tg") = 170,
        py::arg("ng0") = 1.515, py::arg("ng") = 1.515, py::arg("ns") = 1.47,
        py::arg("wvl") = 0.550, py::arg("mag") = 100, py::arg("na") = 1.45,
        py::arg("pixel") = 6.5, py::arg("sf") = 3, py::arg("mode") = 1);

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}

