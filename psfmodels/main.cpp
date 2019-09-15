#include "vectorialPSF.cpp"
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#define cimg_use_tiff
#include "CImg.h"
using namespace cimg_library;

int main(int argc, char *argv[]) {
  //  xp : Source position, 3-element vector [xp yp zp]
  //  nx  : Window size for the psf calculation, in pixels (must be odd).
  //        The origin is located at ((nx+1)/2, (nx+1)/2).
  int nx;
  int nz;
  float range;
  float zdepth; // depth of the point source (positive) from the coverslip
  // p : Parameter structure of system properties, with fields (case sensitive)
  parameters p;

  // Declare the supported options.
  po::options_description desc("Generate vectorial PSF;");
  desc.add_options()("help", "produce help message")(
      "nx", po::value<int>(&nx)->default_value(31),
      "XY size of psf image, in pixels (must be odd)\nThe origin is located at "
      "((nx+1)/2, (nx+1)/2).")("nz", po::value<int>(&nz)->default_value(31),
                               "Z size of psf image, in pixels (must be "
                               "odd)\nThe origin is located at (nz+1)/2.")(
      "range", po::value<float>(&range)->default_value(3.0, "3.0"),
      "Z range of PSF, in microns")(
      "zdepth", po::value<float>(&zdepth)->default_value(0, "0"),
      "position of the point source relative to the coverslip (microns, "
      "positive)")("ti0", po::value<double>(&p.ti0)->default_value(190, "190"),
                   "working distance of the objective (microns)")(
      "ni0", po::value<double>(&p.ni0)->default_value(1.518),
      "immersion medium refractive index, design value")(
      "ni", po::value<double>(&p.ni)->default_value(1.518),
      "immersion medium refractive index, experimental value")(
      "tg0", po::value<double>(&p.tg0)->default_value(170, "170"),
      "coverslip thickness, design value (microns)")(
      "tg", po::value<double>(&p.tg)->default_value(170, "170"),
      "coverslip thickness, experimental value (microns)")(
      "ng0", po::value<double>(&p.ng0)->default_value(1.515, "1.515"),
      "coverslip refractive index, design value")(
      "ng", po::value<double>(&p.ng)->default_value(1.515, "1.515"),
      "coverslip refractive index, experimental value")(
      "ns", po::value<double>(&p.ns)->default_value(1.33, "1.33"),
      "sample refractive index")(
      "lambda", po::value<double>(&p.lambda)->default_value(550),
      "emission wavelength (nm)")(
      "mag", po::value<double>(&p.M)->default_value(100), "magnification")(
      "NA", po::value<double>(&p.NA)->default_value(1.47, "1.47"),
      "Numerical aperture")(
      "pixelSize", po::value<double>(&p.pixelSize)->default_value(6.5, "6.5"),
      "physical size (width) of the camera pixels (microns)")(
      "sf", po::value<int>(&p.sf)->default_value(3),
      "oversampling factor to approximate pixel integration")(
      "mode", po::value<int>(&p.mode)->default_value(1),
      "if 0, returns oversampled PSF");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    cout << desc << "\n";
    return 1;
  }

  p.ti0 = p.ti0 * 1e-6;
  p.tg0 = p.tg0 * 1e-6;
  p.tg = p.tg * 1e-6;
  p.lambda = p.lambda * 1e-9;
  p.ni0_2 = p.ni0 * p.ni0;
  p.ni_2 = p.ni * p.ni;
  p.ng0_2 = p.ng0 * p.ng0;
  p.ng_2 = p.ng * p.ng;
  p.ns_2 = p.ns * p.ns;
  p.k0 = 2 * PI / p.lambda;
  p.NA_2 = p.NA * p.NA;
  p.alpha = asin(p.NA / p.ni);
  p.pixelSize = p.pixelSize * 1e-6;

  zdepth = zdepth * 1e-6;
  double xp[] = {0.0, 0.0, zdepth};
  double lim = range / 2 * 1.e-6;
  std::vector<double> zvec = linspace(-lim + zdepth, lim + zdepth, nz);

  // for (int i = 0; i < zvec.size(); i++) {
  //   std::cout << zvec.at(i) << ' ';
  // }

  VectorialPSF psf = VectorialPSF(xp, &zvec[0], nz, nx, p);
  psf.calculatePSF();
  // psf.calculatePSFdxp();

  int N = nx * nx * nz;
  float out[N];
  std::copy(psf.pixels_, psf.pixels_ + N, out);
  CImg<float> output(nx, nx, nz);
  memcpy(output.data(), out, nx * nx * nz * sizeof(float));

  float dx = p.pixelSize / p.M * 1e6;
  float dz = range / (nz - 1);
  float voxel_size[] = {dx, dx, dz};
  std::string s = "ImageJ=1.50i\n"
                  "spacing=" +
                  std::to_string(dz) +
                  "\n"
                  "unit=micron";
  const char *description = s.c_str();

  output.save_tiff("/Users/talley/Desktop/out.tif", 0, voxel_size, description);

  printf("Done.\n");
  return 0;
}