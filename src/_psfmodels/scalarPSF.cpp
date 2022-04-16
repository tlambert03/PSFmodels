/* scalarPSF.cpp computes the scalar PSF model described by Gibson and Lanni
 * [1]. For more information and implementation details, see [2].
 *
 * [1] F. Aguet et al., Opt. Express 17(8), pp. 6829-6848, 2009
 * [2] F. Aguet, Ph.D Thesis, Swiss Federal Institute of Technology, Lausanne
 * (EPFL), 2009
 *
 * Copyright (C) 2005-2013 Francois Aguet
 *
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * MEX compilation:
 * Mac/Linux: mex -I../../common/mex/include scalarPSF.cpp
 * Windows: mex COMPFLAGS="$COMPFLAGS /MT" -I"..\..\common\mex\include"
 * scalarPSF.cpp
 */

#include "psfmath.h"
#include <string.h>

#define NARGIN 4

using namespace std;

class ScalarPSF
{

public:
  ScalarPSF(const double xp[], const double z[], const int nz, const int nx,
            const parameters p);
  ~ScalarPSF();

  void calculatePSF();
  void calculatePSFdxp();

  double *pixels_;
  double *pixelsDxp_;
  double *pixelsDyp_;
  double *pixelsDzp_;

private:
  double xp_;
  double yp_;
  double zp_;
  const double *z_;
  int nz_;
  int nx_;
  parameters p_;

  int N_;
  double xystep_;

  double **integral_;
  double *R;

  int xymax_;
  int rmax_;
  int npx_;

  static const complex<double> i;
};

const complex<double> ScalarPSF::i = complex<double>(0.0, 1.0);

ScalarPSF::ScalarPSF(const double xp[], const double z[], const int nz,
                     const int nx, const parameters p)
{
  xp_ = xp[0];
  yp_ = xp[1];
  zp_ = xp[2];

  z_ = z;
  nz_ = nz;
  nx_ = nx;
  p_ = p;

  xystep_ = p.dxy;

  xymax_ = ((nx_)*p.sf - 1) / 2; // always fine scale
  if (!p_.mode)
  {
    nx_ *= p_.sf; // oversampling factor
  }

  N_ = nx_ * nx_ * nz_;

  // position in pixels
  xp_ *= p.sf / xystep_;
  yp_ *= p.sf / xystep_;

  int rn = 1 + (int)sqrt(xp_ * xp_ + yp_ * yp_);

  rmax_ = ceil(sqrt(2.0) * xymax_) + rn + 1; // +1 for interpolation, dx, dy
  npx_ = (2 * xymax_ + 1) * (2 * xymax_ + 1);

  pixels_ = new double[N_];
  pixelsDxp_ = new double[N_];
  pixelsDyp_ = new double[N_];
  pixelsDzp_ = new double[N_];

  integral_ = new double *[nz_];
  for (int k = 0; k < nz_; ++k)
  {
    integral_[k] = new double[rmax_];
  }
  // initialize since loops add to these arrays
  memset(pixels_, 0, sizeof(double) * N_);
  memset(pixelsDxp_, 0, sizeof(double) * N_);
  memset(pixelsDyp_, 0, sizeof(double) * N_);
  memset(pixelsDzp_, 0, sizeof(double) * N_);

  // pre-calculate radial coordinates
  R = new double[npx_];
  int idx = 0;
  double xi, yi;
  for (int y = -xymax_; y <= xymax_; ++y)
  {
    for (int x = -xymax_; x <= xymax_; ++x)
    {
      xi = (double)x - xp_;
      yi = (double)y - yp_;
      R[idx] = sqrt(xi * xi + yi * yi);
      ++idx;
    }
  }
}

ScalarPSF::~ScalarPSF()
{
  delete[] R;
  for (int k = 0; k < nz_; ++k)
  {
    delete[] integral_[k];
  }
  delete[] integral_;
  delete[] pixelsDzp_;
  delete[] pixelsDyp_;
  delete[] pixelsDxp_;
  delete[] pixels_;
}

void ScalarPSF::calculatePSF()
{

  double r;
  int n;

  complex<double> sum_I0, expW;

  // constant component of OPD
  double ci = zp_ * (1.0 - p_.ni / p_.ns) +
              p_.ni * (p_.tg0 / p_.ng0 + p_.ti0 / p_.ni0 - p_.tg / p_.ng);

  double theta, sintheta, costheta, ni2sin2theta;
  double bessel_0;

  double A0 = p_.ni_2 * p_.ni_2 / (p_.NA_2 * p_.NA_2);

  // Integration parameters
  double constJ;
  int nSamples;
  double step;

  double w_exp, cst, iconst;
  double ud = 3.0 * p_.sf;

  complex<double> L_th[2];
  for (int k = 0; k < nz_; ++k)
  {

    L_theta(L_th, p_.alpha, p_, ci, z_[k], zp_);
    w_exp = abs(L_th[1]); // missing p.k0, multiply below

    cst = 0.975;
    while (cst >= 0.9)
    {
      L_theta(L_th, cst * p_.alpha, p_, ci, z_[k], zp_);
      if (abs(L_th[1]) > w_exp)
      {
        w_exp = abs(L_th[1]);
      }
      cst -= 0.025;
    }
    w_exp *= p_.k0;

    for (int ri = 0; ri < rmax_; ++ri)
    {
      r = xystep_ / p_.sf * (double)(ri);
      constJ = p_.k0 * r * p_.ni; // samples required for bessel term

      if (w_exp > constJ)
      {
        nSamples = 4 * (int)(1.0 + p_.alpha * w_exp / PI);
      }
      else
      {
        nSamples = 4 * (int)(1.0 + p_.alpha * constJ / PI);
      }
      if (nSamples < 20)
      {
        nSamples = 20;
      }
      step = p_.alpha / (double)nSamples;
      iconst = step / ud;
      iconst *= iconst;

      // Simpson's rule
      sum_I0 = 0.0;
      for (n = 1; n < nSamples / 2; n++)
      {
        theta = 2.0 * n * step;
        sintheta = sin(theta);
        costheta = cos(theta);
        ni2sin2theta = p_.ni_2 * sintheta * sintheta;
        bessel_0 = 2.0 * J0(constJ * sintheta) * sintheta *
                   costheta; // 2.0 factor : Simpson's rule
        expW = exp(i * p_.k0 *
                   ((ci - z_[k]) * p_.ni * costheta +
                    zp_ * sqrt(complex<double>(p_.ns_2 - ni2sin2theta)) +
                    p_.tg * sqrt(complex<double>(p_.ng_2 - ni2sin2theta)) -
                    p_.tg0 * sqrt(complex<double>(p_.ng0_2 - ni2sin2theta)) -
                    p_.ti0 * sqrt(complex<double>(p_.ni0_2 - ni2sin2theta))));
        sum_I0 += expW * bessel_0;
      }
      for (n = 1; n <= nSamples / 2; n++)
      {
        theta = (2.0 * n - 1.0) * step;
        sintheta = sin(theta);
        costheta = cos(theta);
        ni2sin2theta = p_.ni_2 * sintheta * sintheta;
        bessel_0 = 4.0 * J0(constJ * sintheta) * sintheta * costheta;
        expW = exp(i * p_.k0 *
                   ((ci - z_[k]) * p_.ni * costheta +
                    zp_ * sqrt(complex<double>(p_.ns_2 - ni2sin2theta)) +
                    p_.tg * sqrt(complex<double>(p_.ng_2 - ni2sin2theta)) -
                    p_.tg0 * sqrt(complex<double>(p_.ng0_2 - ni2sin2theta)) -
                    p_.ti0 * sqrt(complex<double>(p_.ni0_2 - ni2sin2theta))));
        sum_I0 += expW * bessel_0;
      }
      // theta = alpha;
      bessel_0 = J0(p_.k0 * r * p_.NA) * cos(p_.alpha) * sin(p_.alpha);
      expW = exp(i * p_.k0 *
                 ((ci - z_[k]) * sqrt(complex<double>(p_.ni_2 - p_.NA_2)) +
                  zp_ * sqrt(complex<double>(p_.ns_2 - p_.NA_2)) +
                  p_.tg * sqrt(complex<double>(p_.ng_2 - p_.NA_2)) -
                  p_.tg0 * sqrt(complex<double>(p_.ng0_2 - p_.NA_2)) -
                  p_.ti0 * sqrt(complex<double>(p_.ni0_2 - p_.NA_2))));
      sum_I0 += expW * bessel_0;

      integral_[k][ri] = A0 * abs(sum_I0) * abs(sum_I0) * iconst;
    }
  } // z loop

  int k;
  double dr;

  // Interpolate (linear)
  int r0;
  int index = 0;
  if (p_.mode == 1)
  { // average if sf>1
    div_t divRes;
    for (k = 0; k < nz_; ++k)
    {
      for (int i = 0; i < npx_; ++i)
      {
        r0 = (int)R[i];
        if (r0 + 1 < rmax_)
        {
          dr = R[i] - r0;
          divRes = div(i, 2 * xymax_ + 1);
          index = divRes.rem / p_.sf +
                  (divRes.quot / p_.sf) * nx_; // integer operations!
          pixels_[index + k * nx_ * nx_] +=
              dr * integral_[k][r0 + 1] + (1.0 - dr) * integral_[k][r0];
        } // else '0'
      }
    }
  }
  else
  { // oversample if sf>1
    for (k = 0; k < nz_; ++k)
    {
      for (int i = 0; i < npx_; ++i)
      {
        r0 = (int)R[i];
        if (r0 + 1 < rmax_)
        {
          dr = R[i] - r0;
          pixels_[i + k * npx_] =
              dr * integral_[k][r0 + 1] + (1.0 - dr) * integral_[k][r0];
        } // else '0'
      }
    }
  }
}

void ScalarPSF::calculatePSFdxp()
{

  double r;
  int n;

  double constJ;
  int nSamples;
  double step;

  double theta, sintheta, costheta, ni2sin2theta;
  complex<double> bessel_0, bessel_1, expW, dW, nsroot;
  complex<double> sum_I0, sum_dxI0, sum_dzI0;
  complex<double> tmp;

  // allocate dynamic structures
  double **integralD;
  double **integralDz;
  integralD = new double *[nz_];
  integralDz = new double *[nz_];
  for (int k = 0; k < nz_; k++)
  {
    integralD[k] = new double[rmax_];
    integralDz[k] = new double[rmax_];
  }

  // constant component of OPD
  double ci = zp_ * (1.0 - p_.ni / p_.ns) +
              p_.ni * (p_.tg0 / p_.ng0 + p_.ti0 / p_.ni0 - p_.tg / p_.ng);

  double A0 = p_.ni_2 * p_.ni_2 / (p_.NA_2 * p_.NA_2);

  int ri;
  double ud = 3.0 * p_.sf;

  double w_exp, cst, iconst;

  complex<double> L_th[2];

  for (int k = 0; k < nz_; ++k)
  {

    L_theta(L_th, p_.alpha, p_, ci, z_[k], zp_);
    w_exp = abs(L_th[1]);

    cst = 0.975;
    while (cst >= 0.9)
    {
      L_theta(L_th, cst * p_.alpha, p_, ci, z_[k], zp_);
      if (abs(L_th[1]) > w_exp)
      {
        w_exp = abs(L_th[1]);
      }
      cst -= 0.025;
    }
    w_exp *= p_.k0;

    for (ri = 0; ri < rmax_; ++ri)
    {

      r = xystep_ / p_.sf * (double)(ri);
      constJ = p_.k0 * r * p_.ni;
      if (w_exp > constJ)
      {
        nSamples = 4 * (int)(1.0 + p_.alpha * w_exp / PI);
      }
      else
      {
        nSamples = 4 * (int)(1.0 + p_.alpha * constJ / PI);
      }
      if (nSamples < 20)
      {
        nSamples = 20;
      }
      step = p_.alpha / (double)nSamples;
      iconst = step / ud;
      iconst *= iconst;

      // Simpson's rule
      sum_I0 = 0.0;
      sum_dxI0 = 0.0;
      sum_dzI0 = 0.0;

      for (n = 1; n < nSamples / 2; n++)
      {
        theta = 2.0 * n * step;
        sintheta = sin(theta);
        costheta = cos(theta);
        ni2sin2theta = p_.ni_2 * sintheta * sintheta;
        nsroot = sqrt(complex<double>(p_.ns_2 - ni2sin2theta));

        bessel_0 = 2.0 * J0(constJ * sintheta) * sintheta *
                   costheta; // 2.0 factor : Simpson's rule
        bessel_1 = 2.0 * J1(constJ * sintheta) * sintheta * costheta;

        expW = exp(i * p_.k0 *
                   ((ci - z_[k]) * p_.ni * costheta + zp_ * nsroot +
                    p_.tg * sqrt(complex<double>(p_.ng_2 - ni2sin2theta)) -
                    p_.tg0 * sqrt(complex<double>(p_.ng0_2 - ni2sin2theta)) -
                    p_.ti0 * sqrt(complex<double>(p_.ni0_2 - ni2sin2theta))));
        dW = i * ((1.0 - p_.ni / p_.ns) * p_.ni * costheta + nsroot);

        tmp = expW * bessel_0;
        sum_I0 += tmp;
        tmp *= dW;
        sum_dzI0 += tmp;
        sum_dxI0 += expW * bessel_1 * sintheta;
      }
      for (n = 1; n <= nSamples / 2; n++)
      {
        theta = (2.0 * n - 1.0) * step;
        sintheta = sin(theta);
        costheta = cos(theta);
        ni2sin2theta = p_.ni_2 * sintheta * sintheta;
        nsroot = sqrt(complex<double>(p_.ns_2 - ni2sin2theta));

        bessel_0 = 4.0 * J0(constJ * sintheta) * sintheta *
                   costheta; // 4.0 factor : Simpson's rule
        bessel_1 = 4.0 * J1(constJ * sintheta) * sintheta * costheta;

        expW = exp(i * p_.k0 *
                   ((ci - z_[k]) * p_.ni * costheta + zp_ * nsroot +
                    p_.tg * sqrt(complex<double>(p_.ng_2 - ni2sin2theta)) -
                    p_.tg0 * sqrt(complex<double>(p_.ng0_2 - ni2sin2theta)) -
                    p_.ti0 * sqrt(complex<double>(p_.ni0_2 - ni2sin2theta))));
        dW = i * ((1.0 - p_.ni / p_.ns) * p_.ni * costheta + nsroot);

        tmp = expW * bessel_0;
        sum_I0 += tmp;
        tmp *= dW;
        sum_dzI0 += tmp;
        sum_dxI0 += expW * bessel_1 * sintheta;
      }
      // theta = alpha;
      sintheta = sin(p_.alpha);
      costheta = cos(p_.alpha);
      nsroot = sqrt(complex<double>(p_.ns_2 - p_.NA_2));

      bessel_0 = J0(constJ * sintheta) * sintheta * costheta;
      bessel_1 = J1(constJ * sintheta) * sintheta * costheta;

      expW = exp(i * p_.k0 *
                 ((ci - z_[k]) * p_.ni * costheta + zp_ * nsroot +
                  p_.tg * sqrt(complex<double>(p_.ng_2 - p_.NA_2)) -
                  p_.tg0 * sqrt(complex<double>(p_.ng0_2 - p_.NA_2)) -
                  p_.ti0 * sqrt(complex<double>(p_.ni0_2 - p_.NA_2))));
      dW = i * ((1.0 - p_.ni / p_.ns) * p_.ni * costheta + nsroot);

      tmp = expW * bessel_0;
      sum_I0 += tmp;
      tmp *= dW;
      sum_dzI0 += tmp;
      sum_dxI0 += expW * bessel_1 * sintheta;

      integral_[k][ri] = A0 * abs(sum_I0) * abs(sum_I0) * iconst;
      integralD[k][ri] = p_.k0 * p_.ni * A0 / r * 2.0 *
                         real(conj(sum_I0) * sum_dxI0) *
                         iconst; // multiply with (x-xp)
      integralDz[k][ri] =
          p_.k0 * A0 * 2.0 * real(conj(sum_I0) * sum_dzI0) * iconst;
    }
    integralD[k][0] = 0.0; // overwrite because of singularity
  }                        // z loop

  // Interpolate (linear)
  int r0;
  double dr, rx;
  double xi, yi, tmp2;
  int index = 0;
  int k, x, y;
  if (p_.mode == 1)
  {
    for (k = 0; k < nz_; k++)
    {
      for (y = -xymax_; y <= xymax_; y++)
      {
        for (x = -xymax_; x <= xymax_; x++)
        {

          xi = (double)x - xp_;
          yi = (double)y - yp_;
          rx = sqrt(xi * xi + yi * yi);
          r0 = (int)rx;

          if (r0 + 1 < rmax_)
          {
            dr = rx - r0;
            index = (x + xymax_) / p_.sf + ((y + xymax_) / p_.sf) * nx_ +
                    k * nx_ * nx_;

            pixels_[index] +=
                dr * integral_[k][r0 + 1] + (1.0 - dr) * integral_[k][r0];
            pixelsDzp_[index] +=
                dr * integralDz[k][r0 + 1] + (1.0 - dr) * integralDz[k][r0];

            xi *= xystep_ / p_.sf;
            yi *= xystep_ / p_.sf;

            tmp2 = dr * integralD[k][r0 + 1] + (1.0 - dr) * integralD[k][r0];
            pixelsDxp_[index] += xi * tmp2;
            pixelsDyp_[index] += yi * tmp2;
          } // else '0'
        }
      }
    }
  }
  else
  {
    for (k = 0; k < nz_; k++)
    {
      for (y = -xymax_; y <= xymax_; y++)
      {
        for (x = -xymax_; x <= xymax_; x++)
        {

          xi = (double)x - xp_;
          yi = (double)y - yp_;
          rx = sqrt(xi * xi + yi * yi);
          r0 = (int)rx;

          if (r0 + 1 < rmax_)
          {
            dr = rx - r0;
            pixels_[index] +=
                dr * integral_[k][r0 + 1] + (1.0 - dr) * integral_[k][r0];
            pixelsDzp_[index] +=
                dr * integralDz[k][r0 + 1] + (1.0 - dr) * integralDz[k][r0];

            xi *= xystep_ / p_.sf;
            yi *= xystep_ / p_.sf;

            tmp2 = dr * integralD[k][r0 + 1] + (1.0 - dr) * integralD[k][r0];
            pixelsDxp_[index] += xi * tmp2;
            pixelsDyp_[index] += yi * tmp2;
          } // else '0'
          index++;
        }
      }
    }
  }

  delete[] integralDz;
  delete[] integralD;
}

// compiled with:
// export DYLD_LIBRARY_PATH=/Applications/MATLAB_R2013a.app/bin/maci64 && g++
// -Wall -g -DARRAY_ACCESS_INLINING -I.
// -L/Applications/MATLAB_R2013a.app/bin/maci64 -I../../mex/include/
// -I/Applications/MATLAB_R2013a.app/extern/include scalarPSF.cpp -lmx -lmex
// tested with:
// valgrind --tool=memcheck --leak-check=full --show-reachable=yes ./a.out 2>&1
// | grep scalarPSF

// int main(void) {
//     double xp[] = {0.0, 0.0, 0.0};
//     double z[] = {-100e-9, 0, 100e-9};
//     int nx = 31;
//     int nz = 3;
//
//     parameters p;
//
//     p.ti0 = 1.9e-4;
//     p.ni0 = 1.518;
//     p.ni0_2 = p.ni0*p.ni0;
//     p.ni  = 1.518;
//     p.ni_2 = p.ni*p.ni;
//     p.tg0 = 1.7e-4;
//     p.tg  = 1.7e-4;
//     p.ng0 = 1.515;
//     p.ng0_2 = p.ng0*p.ng0;
//     p.ng  = 1.515;
//     p.ng_2 = p.ng*p.ng;
//     p.ns  = 1.33;
//     p.ns_2 = p.ns*p.ns;
//     p.lambda    = 550e-9;
//     p.k0 = 2*PI/p.lambda;
//     p.M         = 100;
//     p.NA        = 1.45;
//     p.NA_2 = p.NA*p.NA;
//     p.alpha     = asin(p.NA/p.ni);
//     p.pixelSize = 6.45e-6;
//     p.sf = 3;
//     p.mode = 1;
//
//     ScalarPSF psf = ScalarPSF(xp, z, nz, nx, p);
//     psf.calculatePSF();
//     psf.calculatePSFdxp();
//     printf("Done.\n");
// }
