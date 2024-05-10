#include "../VEGAS_interface.hpp"

int main() {
  srand(12);
  double COM = 1000.;
  int nBorn = 3;
  int nUnresolved = 2;
  int reference_index = 2;
  std::vector<int> flavor = {0, 0};
  for(int i = 0 ; i < nBorn; i++) flavor.push_back(1);

  std::vector<Tree<Cluster>> trees = GenTrees(nUnresolved);
  for(int tree_counter = 2; tree_counter < trees.size(); tree_counter++) {
  Tree<Cluster>& tree = trees[tree_counter];
  std::vector<Tree<Cluster>> sectors = GenSectors(flavor, tree, nBorn + 2);
  for(int sec_counter = 1; sec_counter < sectors.size(); sec_counter++) {
  Tree<Cluster> clusterTree = sectors[sec_counter];

  // Test Phase space generation
  std::cout << "\033[1;37m--Test Phase space generation...\033[0m\n" << std::endl;
  PhaseSpace pp = Splitting(nBorn, COM);
  PhaseSpace pp_test = GenMomenta2(pp, clusterTree);
  pp_test.print();
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
  if(!fail) {
    std::cout << "\tMomentum is conserved p = ";
    check.print();
  }

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
  for(int i = 0; i < nUnresolved; i++) {
    std::cout << std::setw(16) << "eta" << i << std::setw(9) << "eta" << i <<" (real)" << std::setw(17) << "Error" << " | ";
  }
  std::cout << std::endl;

  double scale = 1;
  double increment = 0.1;
  bool fail_eta = false;
  while (scale > 1.e-8) {
    scale *= increment;
    std::vector<std::vector<std::vector<double>>> xParFull;
    int level_int = 1;
    std::vector<TreeNode<Cluster>*> level = clusterTree.getLevel(level_int);
    while(level.size() > 0) {
      int unresolved_level = 0;
      for(TreeNode<Cluster>* node : level){
        unresolved_level += node->data.unresolved;

        std::vector<std::vector<double>> xPar;
        for(int c = 0; c < node->data.unresolved; c++) {
          double eta = std::pow(scale, level_int);
          double xi = 0.5;
          double phi = rnd(0., 1.);
          std::vector<double> xPar_c = {eta, xi, phi};
          xPar.push_back(xPar_c);
        }
        xParFull.push_back(xPar);
      }
      level_int++;
      level = clusterTree.getLevel(level_int);
    }
    PhaseSpace ppFull = GenMomenta2(pp, clusterTree, xParFull);
    ppFull.print();
    level_int = 1;
    level = clusterTree.getLevel(level_int);
    while(level.size() > 0) {
      for(TreeNode<Cluster>* node : level){
        Momentum r = node->data.reference_momentum;
        std::cout << "r = "; r.print(); std::cout << std::endl;
        for(Momentum u : node->data.unresolved_momenta) {
          std::cout << "u = "; u.print(); std::cout << std::endl;
          double ratio = r*u/(r.components[0]*u.components[0]*2.)/std::pow(scale, level_int);
          std::cout << std::setw(17) << std::pow(scale, level_int) << std::setprecision(5)
                    << std::setw(17) << r*u/(r.components[0]*u.components[0]*2.) << std::setw(17) << std::abs(1 - ratio) <<  " | ";
          if(std::abs(ratio - 1) > 1.e-5) {
            fail_eta = true;
          }
        }
      }
      level_int++;
      level = clusterTree.getLevel(level_int);
    }
    std::cout << std::endl;
  }


  std::cout << "\n\tLimiting variable xi:" << std::endl;
  for(int i = 0; i < nUnresolved; i++) {
    std::cout << std::setw(16) << "xi" << i << std::setw(9) << "xi" << i <<" (real)" << std::setw(17) << "Error" << " | ";
  }
  std::cout << std::endl;

  scale = 1;
  bool fail_xi = false;
  while (scale > 1.e-8) {
    scale *= increment;
    std::vector<std::vector<std::vector<double>>> xParFull;
    int level_int = 1;
    std::vector<TreeNode<Cluster>*> level = clusterTree.getLevel(level_int);
    while(level.size() > 0) {
      int unresolved_level = 0;
      for(TreeNode<Cluster>* node : level){
        unresolved_level += node->data.unresolved;

        std::vector<std::vector<double>> xPar;
        for(int c = 0; c < node->data.unresolved; c++) {
          double xi = std::pow(scale, level_int);
          double eta = 0.5;
          double phi = rnd(0., 1.);
          std::vector<double> xPar_c = {eta, xi, phi};
          xPar.push_back(xPar_c);
        }
        xParFull.push_back(xPar);
      }
      level_int++;
      level = clusterTree.getLevel(level_int);
    }
    PhaseSpace ppFull = GenMomenta2(pp, clusterTree, xParFull);

    level_int = 1;
    level = clusterTree.getLevel(level_int);
    while(level.size() > 0) {
      for(TreeNode<Cluster>* node : level){
        Momentum r = node->data.reference_momentum;
        for(Momentum u : node->data.unresolved_momenta) {
          double ratio = u.components[0]/(r.components[0] + u.components[0])/std::pow(scale, level_int);
          std::cout << std::setw(17) << std::pow(scale, level_int) << std::setprecision(5)
                    << std::setw(17) << u.components[0]/r.components[0] << std::setw(17) << std::abs(1 - ratio) <<  " | ";
          if(std::abs(ratio - 1) > 1.e-5) {
            fail_xi = true;
          }
          else {
            fail_xi = false;
          }
        }
      }
      level_int++;
      level = clusterTree.getLevel(level_int);
    }
    std::cout << std::endl;
  }

  if(!fail_xi and !fail_eta) {
    std::cout << "\033[1;32m--Infrared limits test successfull!\033[0m\n" << std::endl;
  }
  else {
    std::cout << "\033[1;31m--Infrared limits test failed!\033[0m\n" << std::endl;
  }

  }
  }
}