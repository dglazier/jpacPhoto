// Header file with global phyiscal constants.
// Everything is in GeV unless explicitly stated otherwise.
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _DEBUG_
#define _DEBUG_

#include <iostream>
#include <iomanip>

namespace jpacPhoto
{
  // little function for printing to screen instead of having to copy this line all the time
  template<typename T>
  void debug(T x)
  {
    std::cout << x << std::endl;
  };

  template<typename T, typename F>
  void debug(T x, F y)
  {
    std::cout << std::left << std::setw(15) << x;
    std::cout << std::left << std::setw(15) << y << std::endl;
  };

  template<typename T, typename F, typename G>
  void debug(T x, F y, G z)
  {
    std::cout << std::left << std::setw(15) << x;
    std::cout << std::left << std::setw(15) << y;
    std::cout << std::left << std::setw(15) << z << std::endl;
  };
};

#endif

#ifndef CONSTANT
#define CONSTANT

#include <cmath>
#include <complex>

namespace jpacPhoto
{
  // ---------------------------------------------------------------------------
  const double deg2rad = (M_PI / 180.);
  const double EPS = 1.e-6;
  const double M_ALPHA = 1. / 137.;
  
  const std::complex<double> xr(1., 0.);
  const std::complex<double> xi(0., 1.);
  const std::complex<double> ieps(0., EPS);
  
  // Masses
  const double mPi = 0.13957000; //dg
  const double mK = 0.49367700; //dg
  const double mEta = 0.54785300;
  const double mPro = 0.938272; //dg
  const double mJpsi = 3.0969160;
  const double mPsi2S = 3.686;

  // Masses squared
  const double mPro2 = mPro * mPro;
  const double mJpsi2 = mJpsi * mJpsi;

  // Thresholds
  const double sthPsiPro = (mJpsi + mPro) * (mJpsi + mPro);

  // Decay constants
  const double fJpsi = 0.278;
};
// ---------------------------------------------------------------------------

#endif
