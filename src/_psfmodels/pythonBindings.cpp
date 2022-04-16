#include "scalarPSF.cpp"
#include "vectorialPSF.cpp"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

parameters norm_params(float ti0, float ni0, float ni, float tg0, float tg,
                       float ng0, float ng, float ns, float wvl, float NA,
                       float dxy, int sf, int mode)
{
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
  p.NA = NA;
  p.dxy = dxy * 1e-6;
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

  return p;
}

py::array_t<double> vectorial_psf(py::array_t<double> zv, int nx, double pz,
                                  double ti0, double ni0, double ni, double tg0,
                                  double tg, double ng0, double ng, double ns,
                                  double wvl, double NA, double dxy, int sf = 3,
                                  int mode = 1)
{

  // convert zv microns to meters
  py::buffer_info zvbuf = zv.request();
  double *zvptr = (double *)zvbuf.ptr;
  if (zvbuf.ndim != 1)
    throw std::runtime_error("zv must be a 1-dimensional array");
  for (py::ssize_t i = 0; i < zv.size(); i++)
  {
    zvptr[i] *= 1e-6;
  }

  double xp[] = {0.0, 0.0, pz * 1e-6};

  parameters p =
      norm_params(ti0, ni0, ni, tg0, tg, ng0, ng, ns, wvl, NA, dxy, sf, mode);

  int nz = zv.shape(0);
  VectorialPSF psf = VectorialPSF(xp, zvptr, nz, nx, p);
  psf.calculatePSF();
  return py::array_t<double, py::array::c_style | py::array::forcecast>(
      std::vector<ptrdiff_t>{nz, nx, nx}, &psf.pixels_[0]);
}

py::array_t<double>
vectorial_psf_deriv(py::array_t<double> pixdxp, py::array_t<double> pixdyp,
                    py::array_t<double> pixdzp, py::array_t<double> zv, int nx,
                    double pz, double ti0, double ni0, double ni, double tg0,
                    double tg, double ng0, double ng, double ns, double wvl,
                    double NA, double dxy, int sf = 3, int mode = 1)
{

  // convert zv microns to meters
  py::buffer_info zvbuf = zv.request();
  double *zvptr = (double *)zvbuf.ptr;
  if (zvbuf.ndim != 1)
    throw std::runtime_error("zv must be a 1-dimensional array");
  for (py::ssize_t i = 0; i < zv.size(); i++)
  {
    zvptr[i] *= 1e-6;
  }

  double xp[] = {0.0, 0.0, pz * 1e-6};

  parameters p =
      norm_params(ti0, ni0, ni, tg0, tg, ng0, ng, ns, wvl, NA, dxy, sf, mode);

  int nz = zv.shape(0);
  VectorialPSF psf = VectorialPSF(xp, zvptr, nz, nx, p);
  psf.calculatePSFdxp();

  py::buffer_info pixdxpbuf = pixdxp.request(), pixdypbuf = pixdyp.request(),
                  pixdzpbuf = pixdzp.request();
  double *pixdxpptr = (double *)pixdxpbuf.ptr,
         *pixdypptr = (double *)pixdypbuf.ptr,
         *pixdzpptr = (double *)pixdzpbuf.ptr;

  for (py::ssize_t idx = 0; idx < pixdxpbuf.size; idx++)
  {
    pixdxpptr[idx] = psf.pixelsDxp_[idx];
    pixdypptr[idx] = psf.pixelsDyp_[idx];
    pixdzpptr[idx] = psf.pixelsDzp_[idx];
  }

  return py::array_t<double, py::array::c_style | py::array::forcecast>(
      std::vector<ptrdiff_t>{nz, nx, nx}, &psf.pixels_[0]);
}

py::array_t<double> scalar_psf(py::array_t<double> zv, int nx, double pz,
                               double ti0, double ni0, double ni, double tg0,
                               double tg, double ng0, double ng, double ns,
                               double wvl, double NA, double dxy, int sf = 3,
                               int mode = 1)
{

  // convert zv microns to meters
  py::buffer_info zvbuf = zv.request();
  double *zvptr = (double *)zvbuf.ptr;
  if (zvbuf.ndim != 1)
    throw std::runtime_error("zv must be a 1-dimensional array");
  for (py::ssize_t i = 0; i < zv.size(); i++)
  {
    zvptr[i] *= 1e-6;
  }

  double xp[] = {0.0, 0.0, pz * 1e-6};

  parameters p =
      norm_params(ti0, ni0, ni, tg0, tg, ng0, ng, ns, wvl, NA, dxy, sf, mode);

  int nz = zv.shape(0);
  ScalarPSF psf = ScalarPSF(xp, zvptr, nz, nx, p);
  psf.calculatePSF();
  return py::array_t<double, py::array::c_style | py::array::forcecast>(
      std::vector<ptrdiff_t>{nz, nx, nx}, &psf.pixels_[0]);
}

PYBIND11_MODULE(_psfmodels, m)
{
  m.doc() = "Scalar and Vectorial PSF Models"; // optional module docstring

  m.def("vectorial_psf", &vectorial_psf, R"pbdoc(
    Computes a vectorial microscope point spread function model.

    The model is described in F. Aguet et al., Opt. Express 17(8), pp. 6829-6848, 2009
    For more information and implementation details, see F. Aguet, Ph.D Thesis, Swiss
    Federal Institute of Technology, Lausanne (EPFL), 2009

    C++ code by Francois Aguet, 2009. Python bindings by Talley Lambert, 2019.

    Args:
        zv (np.ndarray): Vector of Z positions at which PSF is calculated (in microns, relative to coverslip)
        nx (int): XY size of output PSF in pixels, must be odd.
        pz (float): point source z position above the coverslip in microns.
        ti0 (float):  working distance of the objective (microns)
        ni0 (float):  immersion medium refractive index, design value
        ni (float):   immersion medium refractive index, experimental value
        tg0 (float):  coverslip thickness, design value (microns)
        tg (float):   coverslip thickness, experimental value (microns)
        ng0 (float):  coverslip refractive index, design value
        ng (float):   coverslip refractive index, experimental value
        ns (float):   sample refractive index
        wvl (float):  emission wavelength (microns)
        NA (float):   numerical aperture
        dxy (float):  pixel size in sample space (microns)
        sf (int):     oversampling factor to approximate pixel integration [default=3]
        mode (int):   if 0, returns oversampled PSF [default=1]

    Returns:
        np.ndarray: PSF with type np.float64 and shape (len(zv), nx, nx)

    )pbdoc",
        py::arg("zv"), py::arg("nx"), py::arg("pz"), py::arg("ti0"),
        py::arg("ni0"), py::arg("ni"), py::arg("tg0"), py::arg("tg"),
        py::arg("ng0"), py::arg("ng"), py::arg("ns"), py::arg("wvl"),
        py::arg("NA"), py::arg("dxy"), py::arg("sf") = 3, py::arg("mode") = 1);

  m.def("vectorial_psf_deriv", &vectorial_psf_deriv, R"pbdoc(
      Computes a vectorial microscope point spread function model, and returns derivatives.
    )pbdoc",
        py::arg("pixdxp"), py::arg("pixdyp"), py::arg("pixdzp"), py::arg("zv"),
        py::arg("nx"), py::arg("pz"), py::arg("ti0"), py::arg("ni0"),
        py::arg("ni"), py::arg("tg0"), py::arg("tg"), py::arg("ng0"),
        py::arg("ng"), py::arg("ns"), py::arg("wvl"), py::arg("NA"),
        py::arg("dxy"), py::arg("sf") = 3, py::arg("mode") = 1);

  m.def("scalar_psf", &scalar_psf, R"pbdoc(
    Computes the scalar PSF model described by Gibson and Lanni

    The model is described in F. Aguet et al., Opt. Express 17(8), pp. 6829-6848, 2009
    For more information and implementation details, see F. Aguet, Ph.D Thesis, Swiss
    Federal Institute of Technology, Lausanne (EPFL), 2009

    C++ code by Francois Aguet, 2009. Python bindings by Talley Lambert, 2019.

    Args:
        zv (np.ndarray): Vector of Z positions at which PSF is calculated (in microns, relative to coverslip)
        nx (int): XY size of output PSF in pixels, must be odd.
        pz (float): point source z position above the coverslip in microns.
        ti0 (float):  working distance of the objective (microns)
        ni0 (float):  immersion medium refractive index, design value
        ni (float):   immersion medium refractive index, experimental value
        tg0 (float):  coverslip thickness, design value (microns)
        tg (float):   coverslip thickness, experimental value (microns)
        ng0 (float):  coverslip refractive index, design value
        ng (float):   coverslip refractive index, experimental value
        ns (float):   sample refractive index
        wvl (float):  emission wavelength (microns)
        NA (float):   numerical aperture
        dxy (float):  pixel size in sample space (microns)
        sf (int):     oversampling factor to approximate pixel integration [default=3]
        mode (int):   if 0, returns oversampled PSF [default=1]

    Returns:
        np.ndarray: PSF with type np.float64 and shape (len(zv), nx, nx)

    )pbdoc",
        py::arg("zv"), py::arg("nx"), py::arg("pz"), py::arg("ti0"),
        py::arg("ni0"), py::arg("ni"), py::arg("tg0"), py::arg("tg"),
        py::arg("ng0"), py::arg("ng"), py::arg("ns"), py::arg("wvl"),
        py::arg("NA"), py::arg("dxy"), py::arg("sf") = 3, py::arg("mode") = 1);
}
