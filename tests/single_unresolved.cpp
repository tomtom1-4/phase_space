#include "../VEGAS_interface.hpp"

int main() {
  srand(12);
  double COM = 1000.;
  int nBorn = 2;
  int nUnresolved = 1;
  int reference_index = 3;
  TreeNode<Cluster>* rootc = new TreeNode<Cluster>(Cluster(reference_index, 1));
  Tree<Cluster> clusterTree;
  clusterTree.setRoot(rootc);

  // Test Phase space generation
  std::cout << "\033[1;37m--Test Phase space generation...\033[0m\n" << std::endl;
  PhaseSpace pp = Splitting(nBorn, COM);
  PhaseSpace pp_test = GenMomenta2(pp, clusterTree);
  clusterTree.print();
  {
  bool fail = false;
  double acc = 1.e-6;
  Momentum check;
  for(int i = 0; i < pp_test.momenta.size(); i++) {
    check = check + pp_test.momenta[i];
  }
  for(int i = 0; i < 4; i++) {
    if(check.components[i] > acc) {
      std::cout << "\tWARNING: Momentum not conserved p = ";
      check.print();
      fail = true;
    }
  }
  std::cout << "\tMomentum is conserved p = ";
  check.print();

  for(int i = 0; i < pp_test.momenta.size(); i++) {
    std::cout << "\tp[" << i << "]^2 = " << pp_test.momenta[i].inv_mass2() << std::endl;
    if (pp_test.momenta[i].inv_mass2() > acc) {
      fail = true;
    }
  }
  if(!fail) {
    std::cout << "\033[1;32m--Phase space generation successfull!\033[0m\n" << std::endl;
  }
  else {
    std::cout << "\033[1;31m--Phase space generation failed!\033[0m\n" << std::endl;
  }
  }
  // Test normalization of Phase space
  int sample_size = 1000000;
  //cubareal x[3*nBorn - 4];
  int spin = -1;
  int neval;
  int fail;
  cubareal integral[1];
  cubareal error[1];
  cubareal prob[1];

  UserData2 data2(COM, nBorn, nUnresolved, &clusterTree);
  std::cout << "\033[1;37m--Perform phase space integration using VEGAS...\033[0m\n" << std::endl;
  Vegas(3*(nBorn + nUnresolved) - 4, 1, *integrand_full2, &data2, 1, 0.001, 0.001, 0, 12, 100, sample_size, 1000, 10000, 1000, 2, "", &spin, &neval, &fail, integral, error, prob);
  std::cout << "\tIntegral = " << integral[0] << " +- " << error[0] << "\t" << error[0]/integral[0]*100 << " %\t" << prob[0] << std::endl;
  PhaseSpace control = RAMBO(nBorn + nUnresolved, COM);
  std::cout << "\tExact result = " << control.weight << std::endl;
  std::cout << "\tRatio = " << integral[0]/control.weight << std::endl;
  std::cout << "\tDeviation = " << (integral[0] - control.weight)/error[0] << std::endl;
  if(std::abs(integral[0] - control.weight)/error[0] < 2.) {
    std::cout << "\033[1;32m--Phase Space Integration successfull!\033[0m\n" << std::endl;
  }
  else {
    std::cout << "\033[1;31m--Phase Space Integration failed!\033[0m\n" << std::endl;
  }

  // Test Infrared limits
  std::cout << "\033[1;37m--Test infrared limits...\033[0m\n" << std::endl;
  std::cout << "\tLimiting variable eta:" << std::endl;
  std::cout << std::setw(17) << "eta" << std::setw(17) << "eta (real)" << std::setw(17) << "Expectation (1)" << std::endl;
  double previous = 0.;
  //pp.print();
  double scale = 1;
  double increment = 0.1;
  bool fail_eta = false;
  while (scale > 1.e-8) {
    scale *= increment;
    std::vector<std::vector<std::vector<double>>> xParFull;
    for(int c = 0; c < nUnresolved; c++) {
      double eta = scale;
      double xi = 0.5;
      double phi = rnd(0., 1.);
      std::vector<std::vector<double>> xPar = {{eta, xi, phi}};
      xParFull.push_back(xPar);
    }
    PhaseSpace ppFull = GenMomenta2(pp, clusterTree, xParFull);
    std::cout << std::setw(17) << scale << std::setw(17) << std::setprecision(8) << ppFull.momenta[reference_index]*ppFull.momenta[reference_index+1]/(ppFull.momenta[reference_index].components[0]*ppFull.momenta[reference_index+1].components[0])/2.
                                               << std::setw(17) << std::setprecision(8) << ppFull.momenta[reference_index]*ppFull.momenta[reference_index+1]/(ppFull.momenta[reference_index].components[0]*ppFull.momenta[reference_index+1].components[0])/2./previous/increment << std::endl;
    previous = ppFull.momenta[reference_index]*ppFull.momenta[reference_index+1]/(ppFull.momenta[reference_index].components[0]*ppFull.momenta[reference_index+1].components[0])/2.;
    if(std::abs(scale - previous)/scale > 1.e-5){
      fail_eta = true;
    }
  }

  std::cout << "\n\tLimiting variable xi:" << std::endl;
  std::cout << std::setw(17) << "xi" << std::setw(17) << "xi (real)" << std::setw(17) << "Expectation (1)" << std::endl;
  scale = 1.;
  bool fail_xi = false;
  while (scale > 1.e-8) {
    scale *= increment;
    std::vector<std::vector<std::vector<double>>> xParFull;
    for(int c = 0; c < nUnresolved; c++) {
      double eta = 0.5;
      double xi = scale;
      double phi = rnd(0., 1.);
      std::vector<std::vector<double>> xPar = {{eta, xi, phi}};
      xParFull.push_back(xPar);
    }
    PhaseSpace ppFull = GenMomenta2(pp, clusterTree, xParFull);
    std::cout << std::setw(17) << scale << std::setw(17) << std::setprecision(8) << ppFull.momenta[reference_index+1].components[0]/(ppFull.momenta[reference_index].components[0])
                                               << std::setw(17) << std::setprecision(8) << ppFull.momenta[reference_index+1].components[0]/(ppFull.momenta[reference_index].components[0])/previous/increment << std::endl;
    previous = ppFull.momenta[reference_index+1].components[0]/(ppFull.momenta[reference_index].components[0]);
    if(std::abs(scale - previous)/scale > 1.e-5){
      fail_xi = true;
    }
    else {
      fail_xi = false;
    }
  }

  if(!fail_xi and !fail_eta) {
    std::cout << "\033[1;32m--Infrared limits test successfull!\033[0m\n" << std::endl;
  }
  else {
    std::cout << "\033[1;31m--Infrared limits test failed!\033[0m\n" << std::endl;
  }

}