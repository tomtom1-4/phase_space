#include "../VEGAS_interface.hpp"

int main() {
  srand(12);
  double COM = 1000.;
  int nBorn = 3;
  int nUnresolved = 2;
  Momentum P(std::vector<double>({COM, 0, 0, 0}));
  TreeNode<Cluster>* root_test = new TreeNode<Cluster>(Cluster(2, 2));
  TreeNode<Cluster>* l1_1 = new TreeNode<Cluster>(Cluster(2, 2));
  TreeNode<Cluster>* l1_2 = new TreeNode<Cluster>(Cluster(2, 2));
  TreeNode<Cluster>* l1_3 = new TreeNode<Cluster>(Cluster(2, 1));
  TreeNode<Cluster>* l2_1_1 = new TreeNode<Cluster>(Cluster(2, 2));
  TreeNode<Cluster>* l2_1_2 = new TreeNode<Cluster>(Cluster(2, 2));


  /*Tree<Cluster> clusterTree_test;
  clusterTree_test.setRoot(root_test);
  clusterTree_test.addChild(root_test, l1_1);
  //clusterTree_test.addChild(root_test, l1_2);
  //clusterTree_test.addChild(root_test, l1_3);
  clusterTree_test.addChild(l1_1, l2_1_1);
  //clusterTree_test.addChild(l1_1, l2_1_2);
  clusterTree_test.print();
  return 0;*/


  {
  int reference_index = 3;
  // One reference momentum, Two unresolved momenta
  std::cout << "\033[1;37m--One reference momentum, Two unresolved momenta:\033[0m\n" << std::endl;
  TreeNode<Cluster>* rootc = new TreeNode<Cluster>(Cluster(reference_index, 2));
  TreeNode<Cluster>* r = new TreeNode<Cluster>(Cluster(reference_index, 0));
  TreeNode<Cluster>* u1 = new TreeNode<Cluster>(Cluster(reference_index, 0));
  TreeNode<Cluster>* u2 = new TreeNode<Cluster>(Cluster(reference_index, 0));

  Tree<Cluster> clusterTree;
  clusterTree.setRoot(rootc);
  clusterTree.addChild(rootc, r);
  clusterTree.addChild(rootc, u1);
  clusterTree.addChild(rootc, u2);
  clusterTree.print();

  // Test Phase space generation
  std::cout << "\033[1;37m--Test Phase space generation...\033[0m\n" << std::endl;
  PhaseSpace pp = Splitting(nBorn, COM);
  PhaseSpace pp_test = GenMomenta2(pp, clusterTree);
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
  std::cout << "\tLimiting variable eta1, eta2:" << std::endl;
  std::cout << std::setw(17) << "scale" << std::setw(17) << "eta1 (real)" << std::setw(17) << "Expectation (1)"
                                       << std::setw(17) << "eta2 (real)" << std::setw(17) << "Expectation (1)" << std::endl;
  double previous1 = 0., previous2 = 0.;
  //pp.print();
  double scale = 1;
  double increment = 0.1;
  bool fail_eta = false;
  while (scale > 1.e-8) {
    scale *= increment;
    std::vector<std::vector<std::vector<double>>> xParFull;
    for(int c = 0; c < 1; c++) {
      double eta1 = scale;
      double xi1 = 0.5;
      double phi1 = rnd(0., 1.);
      double eta2 = scale;
      double xi2 = 0.5;
      double phi2 = rnd(0., 1.);
      std::vector<std::vector<double>> xPar = {{eta1, xi1, phi1}, {eta2, xi2, phi2}};
      xParFull.push_back(xPar);
    }
    PhaseSpace ppFull = GenMomenta2(pp, clusterTree, xParFull);
    std::cout << std::setw(17) << scale << std::setw(17) << std::setprecision(8) << ppFull.momenta[reference_index]*ppFull.momenta[reference_index+1]/(ppFull.momenta[reference_index].components[0]*ppFull.momenta[reference_index+1].components[0])/2.
                                        << std::setw(17) << std::setprecision(8) << ppFull.momenta[reference_index]*ppFull.momenta[reference_index+1]/(ppFull.momenta[reference_index].components[0]*ppFull.momenta[reference_index+1].components[0])/2./previous1/increment
                                        << std::setw(17) << std::setprecision(8) << ppFull.momenta[reference_index]*ppFull.momenta[reference_index+2]/(ppFull.momenta[reference_index].components[0]*ppFull.momenta[reference_index+2].components[0])/2.
                                        << std::setw(17) << std::setprecision(8) << ppFull.momenta[reference_index]*ppFull.momenta[reference_index+2]/(ppFull.momenta[reference_index].components[0]*ppFull.momenta[reference_index+2].components[0])/2./previous2/increment << std::endl;
    previous1 = ppFull.momenta[reference_index]*ppFull.momenta[reference_index+1]/(ppFull.momenta[reference_index].components[0]*ppFull.momenta[reference_index+1].components[0])/2.;
    previous2 = ppFull.momenta[reference_index]*ppFull.momenta[reference_index+2]/(ppFull.momenta[reference_index].components[0]*ppFull.momenta[reference_index+2].components[0])/2.;
    if((std::abs(scale - previous1)/scale > 1.e-5) or (std::abs(scale - previous2)/scale > 1.e-5)){
      fail_eta = true;
    }
  }

  previous1 = 0;
  previous2 = 0;
  std::cout << "\n\tLimiting variable xi1, xi2:" << std::endl;
  std::cout << std::setw(17) << "scale" << std::setw(17) << "xi1 (real)" << std::setw(17) << "Expectation (1)"
                                      << std::setw(17) << "xi2 (real)" << std::setw(17) << "Expectation (1)" << std::endl;
  scale = 1.;
  bool fail_xi = false;
  while (scale > 1.e-8) {
    scale *= increment;
    std::vector<std::vector<std::vector<double>>> xParFull;
    for(int c = 0; c < 1; c++) {
      double eta1 = 0.5;
      double xi1 = scale;
      double phi1 = rnd(0., 1.);
      double eta2 = 0.5;
      double xi2 = scale;
      double phi2 = rnd(0., 1.);
      std::vector<std::vector<double>> xPar = {{eta1, xi1, phi1}, {eta2, xi2, phi2}};
      xParFull.push_back(xPar);
    }
    PhaseSpace ppFull = GenMomenta2(pp, clusterTree, xParFull);
    std::cout << std::setw(17) << scale << std::setw(17) << std::setprecision(8) << ppFull.momenta[reference_index+1].components[0]/(ppFull.momenta[reference_index].components[0])
                                        << std::setw(17) << std::setprecision(8) << ppFull.momenta[reference_index+1].components[0]/(ppFull.momenta[reference_index].components[0])/previous1/increment
                                        << std::setw(17) << std::setprecision(8) << ppFull.momenta[reference_index+2].components[0]/(ppFull.momenta[reference_index].components[0])
                                        << std::setw(17) << std::setprecision(8) << ppFull.momenta[reference_index+2].components[0]/(ppFull.momenta[reference_index].components[0])/previous2/increment << std::endl;
    previous1 = ppFull.momenta[reference_index+1].components[0]/(ppFull.momenta[reference_index].components[0]);
    previous2 = ppFull.momenta[reference_index+2].components[0]/(ppFull.momenta[reference_index].components[0]);
    if((std::abs(scale - previous1)/scale > 1.e-5) or (std::abs(scale - previous2)/scale > 1.e-5)){
      fail_xi = true;
    }
    else {
      fail_xi = false;
    }
  }

  previous1 = 0;
  previous2 = 0;
  std::cout << "\n\tLimiting variable eta1, xi2:" << std::endl;
  std::cout << std::setw(17) << "scale" << std::setw(17) << "eta1 (real)" << std::setw(17) << "Expectation (1)"
                                       << std::setw(17) << "xi2 (real)" << std::setw(17) << "Expectation (1)" << std::endl;
  scale = 1.;
  while (scale > 1.e-8) {
    scale *= increment;
    std::vector<std::vector<std::vector<double>>> xParFull;
    for(int c = 0; c < 1; c++) {
      double eta1 = scale;
      double xi1 = 0.5;
      double phi1 = rnd(0., 1.);
      double eta2 = 0.5;
      double xi2 = scale;
      double phi2 = rnd(0., 1.);
      std::vector<std::vector<double>> xPar = {{eta1, xi1, phi1}, {eta2, xi2, phi2}};
      xParFull.push_back(xPar);
    }
    PhaseSpace ppFull = GenMomenta2(pp, clusterTree, xParFull);
    double u2Max = (P*(pp.momenta[reference_index] - ppFull.momenta[reference_index+1]))/((P - ppFull.momenta[reference_index+1])*ppFull.momenta[reference_index+2])*ppFull.momenta[reference_index+2].components[0];
    std::cout << std::setw(17) << scale << std::setw(17) << std::setprecision(8) << ppFull.momenta[reference_index]*ppFull.momenta[reference_index+1]/(ppFull.momenta[reference_index].components[0]*ppFull.momenta[reference_index+1].components[0])/2.
                                        << std::setw(17) << std::setprecision(8) << ppFull.momenta[reference_index]*ppFull.momenta[reference_index+1]/(ppFull.momenta[reference_index].components[0]*ppFull.momenta[reference_index+1].components[0])/2./previous1/increment
                                        << std::setw(17) << std::setprecision(8) << ppFull.momenta[reference_index+2].components[0]/u2Max
                                        << std::setw(17) << std::setprecision(8) << ppFull.momenta[reference_index+2].components[0]/u2Max/previous2/increment << std::endl;
    previous1 = ppFull.momenta[reference_index]*ppFull.momenta[reference_index+1]/(ppFull.momenta[reference_index].components[0]*ppFull.momenta[reference_index+1].components[0])/2.;
    previous2 = ppFull.momenta[reference_index+2].components[0]/u2Max;
    if((std::abs(scale - previous1)/scale > 1.e-5) or (std::abs(scale - previous2)/scale > 1.e-5)){
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

  {
  int reference_index1 = 2;
  int reference_index2 = 3;
  // Two reference momenta, Two unresolved momenta
  std::cout << "\033[1;37m--Two reference momenta, Two unresolved momenta:\033[0m\n" << std::endl;
  Tree<Cluster> clusterTree2;
  TreeNode<Cluster>* rootc = new TreeNode<Cluster>(Cluster(2, 0));
  TreeNode<Cluster>* node1c = new TreeNode<Cluster>(Cluster(reference_index1, 1));
  TreeNode<Cluster>* node2c = new TreeNode<Cluster>(Cluster(reference_index2, 1));
  TreeNode<Cluster>* r1 = new TreeNode<Cluster>(Cluster(reference_index1, 0));
  TreeNode<Cluster>* u1 = new TreeNode<Cluster>(Cluster(reference_index1+1, 0));
  TreeNode<Cluster>* r2 = new TreeNode<Cluster>(Cluster(reference_index2, 0));
  TreeNode<Cluster>* u2 = new TreeNode<Cluster>(Cluster(reference_index2+1, 0));

  clusterTree2.setRoot(rootc);
  clusterTree2.addChild(rootc, node1c);
  clusterTree2.addChild(rootc, node2c);
  clusterTree2.addChild(node1c, r1);
  clusterTree2.addChild(node2c, r2);
  clusterTree2.addChild(node1c, u1);
  clusterTree2.addChild(node2c, u2);

  Tree<Cluster> clusterTree(clusterTree2);
  clusterTree.print();

  // Test Phase space generation
  std::cout << "\033[1;37m--Test Phase space generation...\033[0m\n" << std::endl;
  PhaseSpace pp = Splitting(nBorn, COM);
  PhaseSpace pp_test = GenMomenta2(pp, clusterTree);
  pp_test.print();

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
  std::cout << "\tLimiting variable eta1, eta2:" << std::endl;
  std::cout << std::setw(17) << "eta1" << std::setw(17) << "eta1 (real)" << std::setw(17) << "Expectation (1)"
            << std::setw(17) << "eta2" << std::setw(17) << "eta2 (real)" << std::setw(17) << "Expectation (1)" << std::endl;
  double previous1 = 0., previous2 = 0.;
  //pp.print();
  double scale = 1;
  double increment = 0.1;
  bool fail_xi = false;
  bool fail_eta = false;
  while (scale > 1.e-5) {
    scale *= increment;
    std::vector<std::vector<std::vector<double>>> xParFull;
    for(int c = 0; c < 1; c++) {
      double eta1 = scale;
      double xi1 = 0.5;
      double phi1 = rnd(0., 1.);
      double eta2 = scale*scale;
      double xi2 = 0.5;
      double phi2 = rnd(0., 1.);
      std::vector<std::vector<double>> xPar1 = {{eta1, xi1, phi1}};
      std::vector<std::vector<double>> xPar2 = {{eta2, xi2, phi2}};
      xParFull.push_back(xPar1);
      xParFull.push_back(xPar2);
    }

    PhaseSpace ppFull = GenMomenta2(pp, clusterTree, xParFull);
    std::cout << std::setw(17) << scale << std::setw(17) << std::setprecision(8) << ppFull.momenta[reference_index1]*ppFull.momenta[reference_index1+1]/(ppFull.momenta[reference_index1].components[0]*ppFull.momenta[reference_index1+1].components[0])/2.
                                        << std::setw(17) << std::setprecision(8) << ppFull.momenta[reference_index1]*ppFull.momenta[reference_index1+1]/(ppFull.momenta[reference_index1].components[0]*ppFull.momenta[reference_index1+1].components[0])/2./previous1/increment
        << std::setw(17) << scale*scale << std::setw(17) << std::setprecision(8) << ppFull.momenta[reference_index2+1]*ppFull.momenta[reference_index2+1+1]/(ppFull.momenta[reference_index2+1].components[0]*ppFull.momenta[reference_index2+1+1].components[0])/2.
                                        << std::setw(17) << std::setprecision(8) << ppFull.momenta[reference_index2+1]*ppFull.momenta[reference_index2+1+1]/(ppFull.momenta[reference_index2+1].components[0]*ppFull.momenta[reference_index2+1+1].components[0])/2./previous2/increment/increment << std::endl;
    previous1 = ppFull.momenta[reference_index1]*ppFull.momenta[reference_index1+1]/(ppFull.momenta[reference_index1].components[0]*ppFull.momenta[reference_index1+1].components[0])/2.;
    previous2 = ppFull.momenta[reference_index2+1]*ppFull.momenta[reference_index2+1+1]/(ppFull.momenta[reference_index2+1].components[0]*ppFull.momenta[reference_index2+1+1].components[0])/2.;
    if((std::abs(scale - previous1)/scale > 1.e-5) or (std::abs(scale*scale - previous2) > 1.e-5)){
      fail_eta = true;
    }
  }

  if(!fail_xi and !fail_eta) {
    std::cout << "\033[1;32m--Infrared limits test successfull!\033[0m\n" << std::endl;
  }
  else {
    std::cout << "\033[1;31m--Infrared limits test failed!\033[0m\n" << std::endl;
  }

  }

}