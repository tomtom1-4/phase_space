#include <iostream>
#include "../src/PhaseSpace.hpp"
#include "../src/Tree.hpp"
#include "../src/Utilities.hpp"


using namespace PSF;
/*
 * This example file illustrates how to generate phase-space points using different algorithms.
 */
int main() {
  double COM = 1000; // Center of Mass Energy
  int nMomenta = 4;
  // Generate random phase space point with RAMBO
  std::cout << "\n################################################################\n" << std::endl;
  std::cout << "RAMBO\n" << std::endl;
  std::cout << "################################################################\n" << std::endl;

  PhaseSpace pp_RAMBO = RAMBO(nMomenta, COM);
  std::cout << "Generated Phase Space Point:" << std::endl;
  pp_RAMBO.print();

  std::cout << "\nCheck Momentum Conservation:" << std::endl;
  pp_RAMBO.check_momentum_conservation();

  std::cout << "\nCheck that momenta are on-shell:" << std::endl;
  pp_RAMBO.check_onshellness();

  std::cout << "\nCorresponding Weight:" << std::endl;
  std::cout << pp_RAMBO.weight << std::endl;

  // Generate random phase space point with Splitting
  std::cout << "\n################################################################\n" << std::endl;
  std::cout << "Splitting (random)\n" << std::endl;
  std::cout << "################################################################\n" << std::endl;

  PhaseSpace pp_Splitting = Splitting(nMomenta, COM);
  std::cout << "Generated Phase Space Point:" << std::endl;
  pp_Splitting.print();

  std::cout << "\nCheck Momentum Conservation:" << std::endl;
  pp_Splitting.check_momentum_conservation();

  std::cout << "\nCheck that momenta are on-shell:" << std::endl;
  pp_Splitting.check_onshellness();


  std::cout << "\nCorresponding Weight:" << std::endl;
  std::cout << pp_Splitting.weight << std::endl;

  // Generate custom phase space point with Splitting
  std::cout << "\n################################################################\n" << std::endl;
  std::cout << "Splitting (custom)\n" << std::endl;
  std::cout << "################################################################\n" << std::endl;

  std::vector<std::vector<double>> x = {{0.5, 0.5} /* invariant masses */,
                                        {0.5, 0.5, 0.5} /* cos */,
                                        {0., 0., 0.} /* phi */};

  pp_Splitting = Splitting(x[0].size() + 2, COM, x);
  std::cout << "Generated Phase Space Point:" << std::endl;
  pp_Splitting.print();

  std::cout << "\nCheck Momentum Conservation:" << std::endl;
  pp_Splitting.check_momentum_conservation();

  std::cout << "\nCheck that momenta are on-shell:" << std::endl;
  pp_Splitting.check_onshellness();


  std::cout << "\nCorresponding Weight:" << std::endl;
  std::cout << pp_Splitting.weight << std::endl;
}