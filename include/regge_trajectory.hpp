// Abstract class for a regge trajectory. Any user-defined class for a specific
// for of the RT can be used in amplitudes (e.g. pomeron_exchange) as long as it
// is derived from this class
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _REGGE_TRAJ_
#define _REGGE_TRAJ_

#include <complex>
#include <string>

class regge_trajectory
{
public:
  // constructor
  regge_trajectory(std::string name = "")
  : parent(name)
  {};

  regge_trajectory(int sig, std::string name = "")
  : signature(sig), parent(name)
  {};

  // copy constructor
  regge_trajectory(const regge_trajectory & old)
  : parent(old.parent), signature(old.signature)
  {};

  // Only need a function to evaluate the trajectory at some s
  virtual std::complex<double> eval(double s) = 0;

  // These parameters define the trajectory
  // name, spin, and mass of the lowest lying resonance on the parent trajectory
  std::string parent;
  int signature;
};


// Basic linear regge_trajectory
class linear_trajectory : public regge_trajectory
{
private:
  // Intercept and slope
  double a0, aprime;

public:
  // Empty constructor
  linear_trajectory(){};

  // Parameterized constructor
  linear_trajectory(int sig, double inter, double slope, std::string name = "")
  : regge_trajectory(sig, name),
    a0(inter), aprime(slope)
  {};

  // copy Constructor
  linear_trajectory(const linear_trajectory & old)
  : regge_trajectory(old),
    a0(old.a0), aprime(old.aprime)
  {};

  // Setting utility
  void set_params(double inter, double slope)
  {
    a0 = inter; aprime = slope;
  };

  std::complex<double> eval(double s)
  {
    return a0 + aprime * s;
  }

  double slope()
  {
    return aprime;
  }
};

#endif
