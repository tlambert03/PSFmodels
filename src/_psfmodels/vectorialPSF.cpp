/* vectorialPSF.cpp computes a vectorial model of the microscope point spread
 * function [1]. For more information and implementation details, see [2].
 *
 * [1] F. Aguet et al., Opt. Express 17(8), pp. 6829-6848, 2009
 * [2] F. Aguet, Ph.D Thesis, Swiss Federal Institute of Technology, Lausanne
 * (EPFL), 2009
 *
 * Copyright (C) 2006-2013 Francois Aguet
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
 * Mac/Linux: mex -I../../mex/include vectorialPSF.cpp
 * Windows: mex COMPFLAGS="$COMPFLAGS /MT" -I"..\..\mex\include"
 * vectorialPSF.cpp
 */

#include "psfmath.h"
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#define NARGIN 4

using namespace std;

template <typename T>
std::vector<T> linspace(T a, T b, size_t N)
{
  T h = (b - a) / static_cast<T>(N - 1);
  std::vector<T> xs(N);
  typename std::vector<T>::iterator x;
  T val;
  for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
    *x = val;
  return xs;
}

class VectorialPSF
{

public:
  VectorialPSF(const double xp[], const double z[], const int nz, const int nx,
               const parameters p);
  ~VectorialPSF();

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

const complex<double> VectorialPSF::i = complex<double>(0.0, 1.0);

VectorialPSF::VectorialPSF(const double xp[], const double z[], const int nz,
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

VectorialPSF::~VectorialPSF()
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

// Intensity PSF for an isotropically emitting point source (average of all
// dipole orientations)
void VectorialPSF::calculatePSF()
{

  double r;
  int n;

  // Integration parameters
  double constJ;
  int nSamples;
  double step;

  double theta, sintheta, costheta, sqrtcostheta, ni2sin2theta;
  complex<double> bessel_0, bessel_1, bessel_2, expW;
  complex<double> ngroot, nsroot;
  complex<double> ts1ts2, tp1tp2;
  complex<double> sum_I0, sum_I1, sum_I2;

  // constant component of OPD
  double ci = zp_ * (1.0 - p_.ni / p_.ns) +
              p_.ni * (p_.tg0 / p_.ng0 + p_.ti0 / p_.ni0 - p_.tg / p_.ng);

  int x, y, index, ri;
  double iconst;
  double ud = 3.0 * p_.sf;

  double w_exp;

  complex<double> L_th[2];
  double cst;

  for (int k = 0; k < nz_; k++)
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

    for (ri = 0; ri < rmax_; ++ri)
    {

      r = xystep_ / p_.sf * (double)(ri);
      constJ = p_.k0 * r * p_.ni; // = w_J;

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

      // Simpson's rule
      sum_I0 = 0.0;
      sum_I1 = 0.0;
      sum_I2 = 0.0;

      for (n = 1; n < nSamples / 2; n++)
      {
        theta = 2.0 * n * step;
        sintheta = sin(theta);
        costheta = cos(theta);
        sqrtcostheta = sqrt(costheta);
        ni2sin2theta = p_.ni_2 * sintheta * sintheta;
        nsroot = sqrt(complex<double>(p_.ns_2 - ni2sin2theta));
        ngroot = sqrt(complex<double>(p_.ng_2 - ni2sin2theta));

        ts1ts2 = 4.0 * p_.ni * costheta * ngroot;
        tp1tp2 = ts1ts2;
        tp1tp2 /= (p_.ng * costheta + p_.ni / p_.ng * ngroot) *
                  (p_.ns / p_.ng * ngroot + p_.ng / p_.ns * nsroot);
        ts1ts2 /= (p_.ni * costheta + ngroot) * (ngroot + nsroot);

        bessel_0 = 2.0 * J0(constJ * sintheta) * sintheta *
                   sqrtcostheta; // 2.0 factor : Simpson's rule
        bessel_1 = 2.0 * J1(constJ * sintheta) * sintheta * sqrtcostheta;
        if (constJ != 0.0)
        {
          bessel_2 = 2.0 * bessel_1 / (constJ * sintheta) - bessel_0;
        }
        else
        {
          bessel_2 = 0.0;
        }
        bessel_0 *= (ts1ts2 + tp1tp2 / p_.ns * nsroot);
        bessel_1 *= (tp1tp2 * p_.ni / p_.ns * sintheta);
        bessel_2 *= (ts1ts2 - tp1tp2 / p_.ns * nsroot);

        expW = exp(i * p_.k0 *
                   ((ci - z_[k]) * p_.ni * costheta + zp_ * nsroot +
                    p_.tg * ngroot -
                    p_.tg0 * sqrt(complex<double>(p_.ng0_2 - ni2sin2theta)) -
                    p_.ti0 * sqrt(complex<double>(p_.ni0_2 - ni2sin2theta))));
        sum_I0 += expW * bessel_0;
        sum_I1 += expW * bessel_1;
        sum_I2 += expW * bessel_2;
      }
      for (n = 1; n <= nSamples / 2; n++)
      {
        theta = (2.0 * n - 1) * step;
        sintheta = sin(theta);
        costheta = cos(theta);
        sqrtcostheta = sqrt(costheta);
        ni2sin2theta = p_.ni_2 * sintheta * sintheta;
        nsroot = sqrt(complex<double>(p_.ns_2 - ni2sin2theta));
        ngroot = sqrt(complex<double>(p_.ng_2 - ni2sin2theta));

        ts1ts2 = 4.0 * p_.ni * costheta * ngroot;
        tp1tp2 = ts1ts2;
        tp1tp2 /= (p_.ng * costheta + p_.ni / p_.ng * ngroot) *
                  (p_.ns / p_.ng * ngroot + p_.ng / p_.ns * nsroot);
        ts1ts2 /= (p_.ni * costheta + ngroot) * (ngroot + nsroot);

        bessel_0 = 4.0 * J0(constJ * sintheta) * sintheta *
                   sqrtcostheta; // 4.0 factor : Simpson's rule
        bessel_1 = 4.0 * J1(constJ * sintheta) * sintheta * sqrtcostheta;
        if (constJ != 0.0)
        {
          bessel_2 = 2.0 * bessel_1 / (constJ * sintheta) - bessel_0;
        }
        else
        {
          bessel_2 = 0.0;
        }
        bessel_0 *= (ts1ts2 + tp1tp2 / p_.ns * nsroot);
        bessel_1 *= (tp1tp2 * p_.ni / p_.ns * sintheta);
        bessel_2 *= (ts1ts2 - tp1tp2 / p_.ns * nsroot);

        expW = exp(i * p_.k0 *
                   ((ci - z_[k]) * p_.ni * costheta + zp_ * nsroot +
                    p_.tg * ngroot -
                    p_.tg0 * sqrt(complex<double>(p_.ng0_2 - ni2sin2theta)) -
                    p_.ti0 * sqrt(complex<double>(p_.ni0_2 - ni2sin2theta))));
        sum_I0 += expW * bessel_0;
        sum_I1 += expW * bessel_1;
        sum_I2 += expW * bessel_2;
      }
      // theta = alpha;
      sintheta = sin(p_.alpha);
      costheta = cos(p_.alpha);
      sqrtcostheta = sqrt(costheta);
      nsroot = sqrt(complex<double>(p_.ns_2 - p_.NA_2));
      ngroot = sqrt(complex<double>(p_.ng_2 - p_.NA_2));

      ts1ts2 = 4.0 * p_.ni * costheta * ngroot;
      tp1tp2 = ts1ts2;
      tp1tp2 /= (p_.ng * costheta + p_.ni / p_.ng * ngroot) *
                (p_.ns / p_.ng * ngroot + p_.ng / p_.ns * nsroot);
      ts1ts2 /= (p_.ni * costheta + ngroot) * (ngroot + nsroot);

      bessel_0 = J0(constJ * sintheta) * sintheta * sqrtcostheta;
      bessel_1 = J1(constJ * sintheta) * sintheta * sqrtcostheta;
      if (constJ != 0.0)
      {
        bessel_2 = 2.0 * bessel_1 / (constJ * sintheta) - bessel_0;
      }
      else
      {
        bessel_2 = 0.0;
      }
      bessel_0 *= (ts1ts2 + tp1tp2 / p_.ns * nsroot);
      bessel_1 *= (tp1tp2 * p_.ni / p_.ns * sintheta);
      bessel_2 *= (ts1ts2 - tp1tp2 / p_.ns * nsroot);

      expW = exp(i * p_.k0 *
                 ((ci - z_[k]) * sqrt(complex<double>(p_.ni_2 - p_.NA_2)) +
                  zp_ * nsroot + p_.tg * ngroot -
                  p_.tg0 * sqrt(complex<double>(p_.ng0_2 - p_.NA_2)) -
                  p_.ti0 * sqrt(complex<double>(p_.ni0_2 - p_.NA_2))));
      sum_I0 += expW * bessel_0;
      sum_I1 += expW * bessel_1;
      sum_I2 += expW * bessel_2;

      sum_I0 = abs(sum_I0);
      sum_I1 = abs(sum_I1);
      sum_I2 = abs(sum_I2);

      integral_[k][ri] =
          8.0 * PI / 3.0 *
          real(sum_I0 * sum_I0 + 2.0 * sum_I1 * sum_I1 + sum_I2 * sum_I2) *
          iconst * iconst;
    }
  } // z loop

  // Interpolate (linear)
  int r0;
  double dr, rx, xi, yi;
  index = 0;
  if (p_.mode == 1)
  {
    for (int k = 0; k < nz_; ++k)
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
          } // else '0'
        }
      }
    }
  }
  else
  {
    for (int k = 0; k < nz_; ++k)
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
          } // else '0'
          index++;
        }
      }
    }
  }
} // psf

// Same PSF calculation as above, but including partial derivatives relative to
// source pos. xp
void VectorialPSF::calculatePSFdxp()
{

  double r;
  int n;

  // Integration parameters
  double constJ;
  int nSamples;
  double step;

  double theta, sintheta, costheta, sqrtcostheta, ni2sin2theta;
  complex<double> bessel_0, bessel_1, bessel_2, bessel_3;
  complex<double> ngroot, nsroot;
  complex<double> ts1ts2, tp1tp2;
  complex<double> sum_I0, sum_I1, sum_I2, sum_dxI0, sum_dxI1, sum_dxI2,
      sum_dzI0, sum_dzI1, sum_dzI2;
  complex<double> t0, t1, t2;
  complex<double> expW, dW, tmp;

  double xystep = p_.dxy;

  // constant component of OPD
  double ci = zp_ * (1.0 - p_.ni / p_.ns) +
              p_.ni * (p_.tg0 / p_.ng0 + p_.ti0 / p_.ni0 - p_.tg / p_.ng);

  // allocate dynamic structures
  double **integralDx;
  double **integralDz;
  integralDx = new double *[nz_];
  integralDz = new double *[nz_];
  for (int k = 0; k < nz_; ++k)
  {
    integralDx[k] = new double[rmax_];
    integralDz[k] = new double[rmax_];
  }

  int x, y, index, ri;
  double iconst;
  double ud = 3.0 * p_.sf;

  double w_exp;

  complex<double> L_th[2];
  double cst;

  for (int k = 0; k < nz_; ++k)
  {

    L_theta(L_th, p_.alpha, p_, ci, z_[k], zp_);
    w_exp = abs(L_th[1]); // missing p.k0 !

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

      r = xystep / p_.sf * (double)(ri);
      constJ = p_.k0 * r * p_.ni; // = w_J;

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

      // Simpson's rule
      sum_I0 = 0.0;
      sum_I1 = 0.0;
      sum_I2 = 0.0;
      sum_dxI0 = 0.0;
      sum_dxI1 = 0.0;
      sum_dxI2 = 0.0;
      sum_dzI0 = 0.0;
      sum_dzI1 = 0.0;
      sum_dzI2 = 0.0;

      for (n = 1; n < nSamples / 2; n++)
      {
        theta = 2.0 * n * step;
        sintheta = sin(theta);
        costheta = cos(theta);
        sqrtcostheta = sqrt(costheta);
        ni2sin2theta = p_.ni_2 * sintheta * sintheta;
        nsroot = sqrt(complex<double>(p_.ns_2 - ni2sin2theta));
        ngroot = sqrt(complex<double>(p_.ng_2 - ni2sin2theta));

        ts1ts2 = 4.0 * p_.ni * costheta * ngroot;
        tp1tp2 = ts1ts2;
        tp1tp2 /= (p_.ng * costheta + p_.ni / p_.ng * ngroot) *
                  (p_.ns / p_.ng * ngroot + p_.ng / p_.ns * nsroot);
        ts1ts2 /= (p_.ni * costheta + ngroot) * (ngroot + nsroot);

        bessel_0 = 2.0 * J0(constJ * sintheta) * sintheta *
                   sqrtcostheta; // 2.0 factor : Simpson's rule
        bessel_1 = 2.0 * J1(constJ * sintheta) * sintheta * sqrtcostheta;
        if (constJ != 0.0)
        {
          bessel_2 = 2.0 * bessel_1 / (constJ * sintheta) - bessel_0;
          bessel_3 = 4.0 * bessel_2 / (constJ * sintheta) - bessel_1;
        }
        else
        {
          bessel_2 = 0.0;
          bessel_3 = 0.0;
        }

        t0 = ts1ts2 + tp1tp2 / p_.ns * nsroot;
        t1 = tp1tp2 * p_.ni / p_.ns * sintheta;
        t2 = ts1ts2 - tp1tp2 / p_.ns * nsroot;

        expW = exp(i * p_.k0 *
                   ((ci - z_[k]) * p_.ni * costheta + zp_ * nsroot +
                    p_.tg * ngroot -
                    p_.tg0 * sqrt(complex<double>(p_.ng0_2 - ni2sin2theta)) -
                    p_.ti0 * sqrt(complex<double>(p_.ni0_2 - ni2sin2theta))));
        dW = i * ((1.0 - p_.ni / p_.ns) * p_.ni * costheta + nsroot);

        tmp = expW * bessel_0 * t0;
        sum_I0 += tmp;
        sum_dzI0 += tmp * dW;
        tmp = expW * bessel_1 * t1;
        sum_I1 += tmp;
        sum_dzI1 += tmp * dW;
        tmp = expW * bessel_2 * t2;
        sum_I2 += tmp;
        sum_dzI2 += tmp * dW;

        sum_dxI0 += expW * bessel_1 * t0 * sintheta;
        sum_dxI1 += expW * (bessel_0 - bessel_2) * t1 * sintheta;
        sum_dxI2 += expW * (bessel_1 - bessel_3) * t2 * sintheta;
      }
      for (n = 1; n <= nSamples / 2; n++)
      {
        theta = (2.0 * n - 1) * step;
        sintheta = sin(theta);
        costheta = cos(theta);
        sqrtcostheta = sqrt(costheta);
        ni2sin2theta = p_.ni_2 * sintheta * sintheta;
        nsroot = sqrt(complex<double>(p_.ns_2 - ni2sin2theta));
        ngroot = sqrt(complex<double>(p_.ng_2 - ni2sin2theta));

        ts1ts2 = 4.0 * p_.ni * costheta * ngroot;
        tp1tp2 = ts1ts2;
        tp1tp2 /= (p_.ng * costheta + p_.ni / p_.ng * ngroot) *
                  (p_.ns / p_.ng * ngroot + p_.ng / p_.ns * nsroot);
        ts1ts2 /= (p_.ni * costheta + ngroot) * (ngroot + nsroot);

        bessel_0 = 4.0 * J0(constJ * sintheta) * sintheta *
                   sqrtcostheta; // 4.0 factor : Simpson's rule
        bessel_1 = 4.0 * J1(constJ * sintheta) * sintheta * sqrtcostheta;
        if (constJ != 0.0)
        {
          bessel_2 = 2.0 * bessel_1 / (constJ * sintheta) - bessel_0;
          bessel_3 = 4.0 * bessel_2 / (constJ * sintheta) - bessel_1;
        }
        else
        {
          bessel_2 = 0.0;
          bessel_3 = 0.0;
        }
        t0 = ts1ts2 + tp1tp2 / p_.ns * nsroot;
        t1 = tp1tp2 * p_.ni / p_.ns * sintheta;
        t2 = ts1ts2 - tp1tp2 / p_.ns * nsroot;

        expW = exp(i * p_.k0 *
                   ((ci - z_[k]) * p_.ni * costheta + zp_ * nsroot +
                    p_.tg * ngroot -
                    p_.tg0 * sqrt(complex<double>(p_.ng0_2 - ni2sin2theta)) -
                    p_.ti0 * sqrt(complex<double>(p_.ni0_2 - ni2sin2theta))));
        dW = i * ((1.0 - p_.ni / p_.ns) * p_.ni * costheta + nsroot);

        tmp = expW * bessel_0 * t0;
        sum_I0 += tmp;
        sum_dzI0 += tmp * dW;
        tmp = expW * bessel_1 * t1;
        sum_I1 += tmp;
        sum_dzI1 += tmp * dW;
        tmp = expW * bessel_2 * t2;
        sum_I2 += tmp;
        sum_dzI2 += tmp * dW;

        sum_dxI0 += expW * bessel_1 * t0 * sintheta;
        sum_dxI1 += expW * (bessel_0 - bessel_2) * t1 * sintheta;
        sum_dxI2 += expW * (bessel_1 - bessel_3) * t2 * sintheta;
      }
      // theta = alpha;
      sintheta = sin(p_.alpha);
      costheta = cos(p_.alpha);
      sqrtcostheta = sqrt(costheta);
      nsroot = sqrt(complex<double>(p_.ns_2 - p_.NA_2));
      ngroot = sqrt(complex<double>(p_.ng_2 - p_.NA_2));

      ts1ts2 = 4.0 * p_.ni * costheta * ngroot;
      tp1tp2 = ts1ts2;
      tp1tp2 /= (p_.ng * costheta + p_.ni / p_.ng * ngroot) *
                (p_.ns / p_.ng * ngroot + p_.ng / p_.ns * nsroot);
      ts1ts2 /= (p_.ni * costheta + ngroot) * (ngroot + nsroot);

      bessel_0 = J0(constJ * sintheta) * sintheta * sqrtcostheta;
      bessel_1 = J1(constJ * sintheta) * sintheta * sqrtcostheta;
      if (constJ != 0.0)
      {
        bessel_2 = 2.0 * bessel_1 / (constJ * sintheta) - bessel_0;
        bessel_3 = 4.0 * bessel_2 / (constJ * sintheta) - bessel_1;
      }
      else
      {
        bessel_2 = 0.0;
        bessel_3 = 0.0;
      }
      t0 = ts1ts2 + tp1tp2 / p_.ns * nsroot;
      t1 = tp1tp2 * p_.ni / p_.ns * sintheta;
      t2 = ts1ts2 - tp1tp2 / p_.ns * nsroot;

      expW = exp(i * p_.k0 *
                 ((ci - z_[k]) * sqrt(complex<double>(p_.ni_2 - p_.NA_2)) +
                  zp_ * nsroot + p_.tg * ngroot -
                  p_.tg0 * sqrt(complex<double>(p_.ng0_2 - p_.NA_2)) -
                  p_.ti0 * sqrt(complex<double>(p_.ni0_2 - p_.NA_2))));
      dW = i * ((1.0 - p_.ni / p_.ns) * p_.ni * costheta + nsroot);

      tmp = expW * bessel_0 * t0;
      sum_I0 += tmp;
      sum_dzI0 += tmp * dW;
      tmp = expW * bessel_1 * t1;
      sum_I1 += tmp;
      sum_dzI1 += tmp * dW;
      tmp = expW * bessel_2 * t2;
      sum_I2 += tmp;
      sum_dzI2 += tmp * dW;

      sum_dxI0 += expW * bessel_1 * t0 * sintheta;
      sum_dxI1 += expW * (bessel_0 - bessel_2) * t1 * sintheta;
      sum_dxI2 += expW * (bessel_1 - bessel_3) * t2 * sintheta;

      if (ri > 0)
      {
        integral_[k][ri] =
            8.0 * PI / 3.0 *
            (abs(sum_I0) * abs(sum_I0) + 2.0 * abs(sum_I1) * abs(sum_I1) +
             abs(sum_I2) * abs(sum_I2)) *
            iconst * iconst;
        integralDx[k][ri] =
            16.0 * PI / 3.0 * p_.k0 * p_.ni *
            real(-sum_dxI0 * conj(sum_I0) + sum_dxI1 * conj(sum_I1) +
                 sum_dxI2 * conj(sum_I2) / 2.0) /
            r * iconst * iconst;
        integralDz[k][ri] =
            16.0 * PI / 3.0 * p_.k0 *
            real(conj(sum_dzI0) * sum_I0 + 2.0 * conj(sum_dzI1) * sum_I1 +
                 conj(sum_dzI2) * sum_I2) *
            iconst * iconst;
      }
      else
      {
        integral_[k][0] =
            8.0 * PI / 3.0 * (abs(sum_I0) * abs(sum_I0)) * iconst * iconst;
        integralDx[k][0] = 0.0;
        integralDz[k][0] = 16.0 * PI / 3.0 * p_.k0 *
                           real(sum_I0 * conj(sum_dzI0)) * iconst * iconst;
      }
    }
  } // z loop

  // Interpolate (linear)
  int r0;
  double dr, rx, xi, yi, xd, yd;
  index = 0;
  if (p_.mode == 1)
  {
    for (int k = 0; k < nz_; ++k)
    {
      for (y = -xymax_; y <= xymax_; y++)
      {
        for (x = -xymax_; x <= xymax_; x++)
        {
          xi = (double)x - xp_;
          yi = (double)y - yp_;
          xd = xp_ - x * xystep / p_.sf;
          yd = yp_ - y * xystep / p_.sf;
          rx = sqrt(xi * xi + yi * yi);
          r0 = (int)rx;
          if (r0 + 1 < rmax_)
          {
            dr = rx - r0;
            index = (x + xymax_) / p_.sf + ((y + xymax_) / p_.sf) * nx_ +
                    k * nx_ * nx_;
            pixels_[index] +=
                dr * integral_[k][r0 + 1] + (1.0 - dr) * integral_[k][r0];
            pixelsDxp_[index] += xd * (dr * integralDx[k][r0 + 1] +
                                       (1.0 - dr) * integralDx[k][r0]);
            pixelsDyp_[index] += yd * (dr * integralDx[k][r0 + 1] +
                                       (1.0 - dr) * integralDx[k][r0]);
            pixelsDzp_[index] +=
                dr * integralDz[k][r0 + 1] + (1.0 - dr) * integralDz[k][r0];
          } // else '0'
        }
      }
    }
  }
  else
  {
    for (int k = 0; k < nz_; ++k)
    {
      for (y = -xymax_; y <= xymax_; y++)
      {
        for (x = -xymax_; x <= xymax_; x++)
        {
          xi = (double)x - xp_;
          yi = (double)y - yp_;
          xd = xp_ - x * xystep_ / p_.sf;
          yd = yp_ - y * xystep_ / p_.sf;
          rx = sqrt(xi * xi + yi * yi);
          r0 = (int)rx;
          if (r0 + 1 < rmax_)
          {
            dr = rx - r0;
            pixels_[index] +=
                dr * integral_[k][r0 + 1] + (1.0 - dr) * integral_[k][r0];
            pixelsDxp_[index] += xd * (dr * integralDx[k][r0 + 1] +
                                       (1.0 - dr) * integralDx[k][r0]);
            pixelsDyp_[index] += yd * (dr * integralDx[k][r0 + 1] +
                                       (1.0 - dr) * integralDx[k][r0]);
            pixelsDzp_[index] +=
                dr * integralDz[k][r0 + 1] + (1.0 - dr) * integralDz[k][r0];
          } // else '0'
          index++;
        }
      }
    }
  }
  // free dynamic structures
  for (int k = 0; k < nz_; ++k)
  {
    delete[] integralDx[k];
    delete[] integralDz[k];
  }
  delete[] integralDx;
  delete[] integralDz;
} // psf_dx

// [h, dxp, dyp, dzp] = vectorialPSF(xp, z, nx, p) computes a vectorial
// microscope point spread function model.
//  The partial derivatives of the model relative to the source position xp are
//  also calculated. The model is described in [1]. For more information and
//  implementation details, see [2].
//
//    INPUTS:
//    xp    : Source position, 3-element vector [xp yp zp]
//    z     : Vector of z-plane positions
//    nx    : Window size for the psf calculation, in pixels (must be odd).
//            The origin is located at ((nx+1)/2, (nx+1)/2).
//    p     : Parameter structure of system properties, with fields (case
//    sensitive)
//             ti0       : working distance of the objective
//             ni0       : immersion medium refractive index, design value
//             ni        : immersion medium refractive index, experimental value
//             tg0       : coverslip thickness, design value
//             tg        : coverslip thickness, experimental value
//             ng0       : coverslip refractive index, design value
//             ng        : coverslip refractive index, experimental value
//             ns        : sample refractive index
//             lambda    : emission wavelength
//             M         : magnification
//             NA        : numerical aperture
//             pixelSize : physical size (width) of the camera pixels
//             f         : (optional, default: 3) oversampling factor to
//             approximate pixel integration mode      : (optional, default: 1)
//             if 0, returns oversampled PSF
//
//    All spatial units are in object space, in [m].

// void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// {

//     // Input checks
//     // if (nrhs!=NARGIN) mexErrMsgTxt("There must be 4 input arguments: xp,
//     z, w, p");
//     // if ( !mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) ||
//     !mxIsDouble(prhs[2]) ) mexErrMsgTxt("'xp' and 'z' must be double arrays.
//     Window size 'w' must be an integer.");
//     // if ( !mxIsDouble(prhs[2]) ) mexErrMsgTxt("Input 'z' must be a double
//     array.");
//     // if ( !mxIsStruct(prhs[3]) ) mexErrMsgTxt("Input 'p' must be a
//     parameter structure.");

//     if (mxGetNumberOfElements(prhs[0])!=3) mexErrMsgTxt("Input 'xp' must be a
//     3-element vector."); double* xp = mxGetPr(prhs[0]);

//     int nz = (int)mxGetNumberOfElements(prhs[1]);
//     double* z = mxGetPr(prhs[1]);

//     int nx = (int)mxGetScalar(prhs[2]);
//     if (nx%2!=1) {
//         mexErrMsgTxt("Windows size must be an odd integer");
//     }

//     int np = mxGetNumberOfFields(prhs[3]);
//     if (np < 12) mexErrMsgTxt("Incorrect parameter vector");
//     parameters p;
//     parseParameterStruct(p, prhs[3]);

//     VectorialPSF psf = VectorialPSF(xp, z, nz, nx, p);
//     if (nlhs==1) {
//         psf.calculatePSF();
//     } else if (nlhs>1) {
//         psf.calculatePSFdxp();
//     }

//     int ndim = 3;
//     if (p.mode==0) {
//         nx *= p.sf;
//     }
//     const mwSize dims[3] = {nx, nx, nz};

//     int N = nx*nx*nz;
//     // copy PSF data
//     plhs[0] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
//     memcpy(mxGetPr(plhs[0]), psf.pixels_, N*sizeof(double));

//     // copy derivatives
//     if (nlhs>1) {
//         plhs[1] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
//         memcpy(mxGetPr(plhs[1]), psf.pixelsDxp_, N*sizeof(double));
//     }
//     if (nlhs>2) {
//         plhs[2] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
//         memcpy(mxGetPr(plhs[2]), psf.pixelsDyp_, N*sizeof(double));
//     }
//     if (nlhs>3) {
//         plhs[3] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
//         memcpy(mxGetPr(plhs[3]), psf.pixelsDzp_, N*sizeof(double));
//     }
// }

// compiled with:
// export DYLD_LIBRARY_PATH=/Applications/MATLAB_R2013a.app/bin/maci64 && g++
// -Wall -g -DARRAY_ACCESS_INLINING -I.
// -L/Applications/MATLAB_R2013a.app/bin/maci64 -I../../mex/include/
// -I/Applications/MATLAB_R2013a.app/extern/include vectorialPSF.cpp -lmx -lmex
// tested with:
// valgrind --tool=memcheck --leak-check=full --show-reachable=yes ./a.out 2>&1
// | grep vectorialPSF
