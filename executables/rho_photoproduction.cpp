// ---------------------------------------------------------------------------
// Photoproduction of rho by a  charged pion exchange
// Reproduces the results from [1]
//
// Author:       Derek Glazier (2020)
// ---------------------------------------------------------------------------
// COMMAND LINE OPTIONS:
// -c double          # Change CM angle in degree (default: 0)
// -n int             # Number of points to plot (default: 25)
// -m double          # Maximum CM angle to plot (default: 10 GeV)
// -diff              # Plot differential xsection (default: false)
// -y "[y1:y2]"       # Custom y bounds in output plot
// -lab               # Display E_lab in the x-axis (default: false)
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "regge_trajectory.hpp"
#include "amplitudes/pseudoscalar_exchange.hpp"
#include "amplitudes/pomeron_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

#include <TRandom3.h>
#include <TH1F.h>
#include <TF2.h>
#include <TBenchmark.h>
#include <Math/MinimizerOptions.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>

#include <cstring>
#include <iostream>

using namespace jpacPhoto;

int main( int argc, char** argv )
{

  // ---------------------------------------------------------------------------
  // COMMAND LINE OPTIONS
  // ---------------------------------------------------------------------------

  double theta = 0.;
  double max = 50;
  double y[2]; bool custom_y = false;
  int N = 25;
  std::string xlabel = "W   [GeV]"; bool LAB = false;
  std::string ylabel = "#sigma(#gamma N #rightarrow #rho N)   [ub]";
  std::string filename = "rho_photoproduction_sum.pdf";
  bool INTEG = true;

  // ---------------------------------------------------------------------------
  // Parse command line arguments
  for (int i = 0; i < argc; i++)
  {
    if (std::strcmp(argv[i],"-f")==0) filename = argv[i+1];
    if (std::strcmp(argv[i],"-c")==0) theta = atof(argv[i+1]);
    if (std::strcmp(argv[i],"-m")==0) max = atof(argv[i+1]);
    if (std::strcmp(argv[i],"-n")==0) N = atoi(argv[i+1]);
    if (std::strcmp(argv[i],"-y")==0)
    {
      custom_y = true;
      y_range(argv[i+1], y);
    }
    if (std::strcmp(argv[i],"-diff")==0)
    {
      INTEG = false;
      ylabel = "d#sigma/dt  [#mub GeV^{-2}]";
    }
    if (std::strcmp(argv[i],"-lab")==0)
    {
      LAB = true;
      xlabel = "E_{#gamma}   [GeV]";
    }
  }

  // ---------------------------------------------------------------------------
  // AMPLITUDES
  // ---------------------------------------------------------------------------


  // ---------------------------------------------------------------------------
  // rho

  // Kinematics for 780 MeV vector
  double rhoMass=0.77549000;
  double rhoWidth=0.15100000;
  reaction_kinematics * ptrRho = new reaction_kinematics(rhoMass, "rho0");

  // Amplitudes
  // For regge amps need pion trajectory
  linear_trajectory alpha(1, -0.7*mPi*mPi, 0.7, "pionic trajectory");
  pseudoscalar_exchange rho_pi(ptrRho, mPi, "#pi");
  pseudoscalar_exchange rho_regge(ptrRho, &alpha, "Reggeon exchange");
  // Couplings for rho  width
  rho_pi.set_params({5, sqrt(rhoWidth*1000.*M_PI*14.4)});
  rho_regge.set_params({4, sqrt(rhoWidth*1000.*M_PI*14.4)});

  // Best fit values from [1] from high energy
  linear_trajectory alpha2016(+1, 1.1, 0.11, "pomeron");
  pomeron_exchange rho_pomeron(ptrRho, &alpha2016, false, "pomeron");
  rho_pomeron.set_params({30,1});
  //rho_pomeron.set_params({30,sqrt(rhoWidth*1000.*M_PI*14.4)});

  //rho_pi.set_params({0.67 * 3.90, sqrt(rhoWidth*1000.*M_PI*14.4)});
  //rho_regge.set_params({0.67 * 3.90, sqrt(rhoWidth*1000.*M_PI*14.4)});

  amplitude_sum sum(ptrRho, {&rho_pi, &rho_pomeron}, "Sum");
  // Add to a vector to plot them both
  std::vector<amplitude*> amps;
  //  amps.push_back(&rho_pi);
  //amps.push_back(&rho_pomeron);
  amps.push_back(&sum);
 
  // ---------------------------------------------------------------------------
  // You shouldnt need to change anything below this line
  // ---------------------------------------------------------------------------

  // Plotter object
  jpacGraph1D* plotter = new jpacGraph1D();
  
  // ---------------------------------------------------------------------------
  // Print the desired observable for each amplitude
  for (int n = 0; n < amps.size(); n++)
  {
    std::cout << std::endl << "Printing amplitude: " << amps[n]->identifier << "\n";

    double low;
    (LAB == true) ? (low = E_lab(amps[n]->kinematics->Wth) + EPS)
                  : (low = amps[n]->kinematics->Wth + EPS);
    auto F = [&](double x)
    {
      double s;
      (LAB == false) ? (s = x*x) : (s = W_cm(x) * W_cm(x));

      if (INTEG == false)
      {
        double t = amps[n]->kinematics->t_man(s, theta * deg2rad);
        return amps[n]->differential_xsection(s, t)/1000;
      }
      else
      {
        return amps[n]->integrated_xsection(s)/1000;
      }
    };


    std::array<std::vector<double>, 2> x_fx = vec_fill(N, F, low, max, true);
    plotter->AddEntry(x_fx[0], x_fx[1], amps[n]->identifier);
  }

  double low;
  (LAB == true) ? (low = E_lab(ptrRho->Wth) + EPS)
                : (low = ptrRho->Wth + EPS);
  plotter->SetXaxis(xlabel, low, max);

  // To change the range of the Y-axis or the position of the Legend change the arguments here
  (custom_y == true) ? (plotter->SetYaxis(ylabel, y[0], y[1])) : (plotter->SetYaxis(ylabel));

  // Position of the legend
  plotter->SetLegend(0.2, 0.65);

  // Output to file
  plotter->Plot(filename);

  delete ptrRho, plotter;
  return 1.;
}
