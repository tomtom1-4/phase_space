#include <iostream>
#include "../src/PhaseSpace.hpp"
#include "../src/Tree.hpp"
#include "../src/Utilities.hpp"

#ifdef USE_CUBA
#include "/home/tom/Documents/software/software/Cuba-4.2.2/cuba.h"
#include "../VEGAS_interface.hpp"
#endif

using namespace PSF;
/*
 * This example file shows how to use the phase-space generation algorithms to perform general phase space integrations.
 * Here, we focus on the integration over the identity as it allows for a straight forward validation.
 */
int main() {
  double COM = 1000; // Center of Mass Energy
  int nMomenta = 4; // Number of final state momenta to integrate over

  // Since RAMBO has perfect efficiency, it has constant weight and therefore the weight directly gives the total phase-space volume.
  PhaseSpace pp_RAMBO = RAMBO(nMomenta, COM);
  std::cout << "Phase-space volume:" << std::endl;
  double control = pp_RAMBO.weight;
  std::cout << control << std::endl;

  // Now we try to compute the phase space volume using different phase-space generation algorithms
  std::cout << "\n################################################################\n" << std::endl;
  std::cout << "Splitting\n" << std::endl;
  std::cout << "################################################################\n" << std::endl;
  int samples = 100000;
  double var = 0.;
  double err = 0.;
  double result = 0.;
  for(int event_counter = 0; event_counter < samples; event_counter++) {
    PhaseSpace pp = Splitting(nMomenta, COM);
    double integrand = 1.;

    result = result*event_counter/(event_counter + 1.) + pp.weight*integrand/(event_counter + 1.);
    var = var*event_counter/(event_counter + 1.) + std::pow(pp.weight*integrand - result, 2)/(event_counter + 1.);

    err = std::sqrt(var)/std::sqrt(event_counter);


    if((event_counter + 1) % (samples/10) == 0) {
      std::cout << std::setw(10) << std::setprecision(5) << event_counter + 1 << "/" << samples
                  << std::setw(15) << result
                  << " +- " << std::setw(6) << err
                  << std::setw(2) << "(" << std::setw(7) << err/result*100 << "%)"
                  << "    deviation: " << std::setw(10) << (result - control)/err
                  << "    ratio: " << std::setw(10) << result/control
                  << std::endl;
    }
  }

  #ifdef USE_CUBA
  std::cout << "\nPerform integration using the VEGAS algorithm with CUBA" << std::endl;
  int spin = -1;
  int neval;
  int fail;
  cubareal integral[1];
  cubareal error[1];
  cubareal prob[1];
  Vegas(3*nMomenta - 4, 1, *integrand_born, &COM, 1, 0.0001, 0.0001, 8, 12, 100, samples, 1000, 10000, 1000, 1, "", &spin, &neval, &fail, integral, error, prob);
  std::cout << "integral = " << integral[0] << " +- " << error[0] << " (" << error[0]/integral[0]*100 << "%)\t probability that error is accurate: " << prob[0] << std::endl;
  #endif
}