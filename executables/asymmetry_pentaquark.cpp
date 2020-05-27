// ---------------------------------------------------------------------------
// Predicted sensativity to LHCb pentaquarks in beam and parity asymmetries at
// GlueX at JLab
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
// References:
// [1] 10.1103/PhysRevD.100.034019
// [2] 10.1103/PhysRevLett.115.072001
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "amplitudes/baryon_resonance.hpp"
#include "amplitudes/pomeron_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"

#include "jpacGraph1D.hpp"

#include <cstring>
#include <iostream>
#include <iomanip>

int main( int argc, char** argv )
{
  double egam = 9.;
  std::string filename = "5q_beam_asymmetry.pdf";
  for (int i = 0; i < argc; i++)
  {
    if (std::strcmp(argv[i],"-e")==0) egam = atof(argv[i+1]);
    if (std::strcmp(argv[i],"-f")==0) filename = argv[i+1];
  }

  // Set up Kinematics
  reaction_kinematics * ptr = new reaction_kinematics(mJpsi, "jpsi");

  // Incoherent sum of the s and t channels
  amplitude_sum sum(ptr, "Sum");

  // ---------------------------------------------------------------------------
  // S - CHANNEL

  // Two different pentaquarks
  // masses and widths from 2015 LHCb paper [2]
  baryon_resonance P_c1(ptr, 3, -1, 4.45, 0.040, "P_{c}(4450)");
  baryon_resonance P_c2(ptr, 5, +1, 4.38, 0.205, "P_{c}(4380)");

  // 1% branching fraction and equal photocouplings for both
  std::vector<double> params = {0.01, .7071};
  P_c1.set_params(params);
  P_c2.set_params(params);

  // Add them to the sum
  sum.add_amplitude(&P_c1);
  sum.add_amplitude(&P_c2);

  // ---------------------------------------------------------------------------
  // T - CHANNEL

  // Set up pomeron trajectory
  linear_trajectory alpha(+1, 0.941, 0.364);

  // Create amplitude with kinematics and trajectory
  pomeron_exchange background(ptr, &alpha, "Background");

  // normalization and t-slope
  std::vector<double> back_params = {0.367, 0.12};
  background.set_params(back_params);

  // Add to the sum
  sum.add_amplitude(&background);

  std::vector<amplitude*> amps = {&sum, &background, &P_c1, &P_c2};

  int N = 200; // how many points to plot

// ---------------------------------------------------------------------------
// You shouldnt need to change anything below this line
// ---------------------------------------------------------------------------

// Plotter objects
jpacGraph1D* plotter = new jpacGraph1D();

// ---------------------------------------------------------------------------
// scan over energy
double s = 2.*mPro* egam + mPro2;

for (int n = 0; n < amps.size(); n++)
{
  std::vector<double> theta, sigma;
  for (int i = 1; i <= N; i++)
  {
    double theta_i = double(i) * 180. / N;
    theta.push_back(theta_i);

    double sigma_i = amps[n]->beam_asymmetry(s, cos(theta_i * deg2rad));
    sigma.push_back(sigma_i);
  }

  plotter->AddEntry(theta, sigma, amps[n]->identifier);
}

std::ostringstream streamObj;
streamObj << std::setprecision(3);
streamObj << egam;

plotter->SetLegend(0.18, 0.7, "E_{#gamma} = " + streamObj.str() + " GeV");
plotter->SetXaxis("#theta  (GeV)", 0., 180.);
plotter->SetYaxis("#Sigma_{4#pi}", -0.4, 1.);

plotter->Plot(filename.c_str());

return 1.;
};