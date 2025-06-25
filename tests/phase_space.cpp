#include "../VEGAS_interface.hpp"

using namespace PSF;

int main() {
  srand(12);
  double COM = 1000.;
  int nBorn = 3;
  int nUnresolved = 2;
  std::vector<bool> flavor = {0, 0};
  for(int i = 0 ; i < nBorn; i++) flavor.push_back(1);

  std::vector<Tree<Cluster>> trees = GenTrees(nUnresolved);
  for(int tree_counter = 0; tree_counter < trees.size(); tree_counter++) {
  Tree<Cluster>& tree = trees[tree_counter];
  if(tree.getRoot()->children.size() > 1) continue; // Remove multiple reference, as they are trivial
  std::vector<Tree<Cluster>> sectors = GenSectors(flavor, tree, nBorn + 2);
  //for(int sec_counter = 0; sec_counter < sectors.size(); sec_counter++) {
  for(int sec_counter = 0; sec_counter < 1; sec_counter++) {
  Tree<Cluster> clusterTree = sectors[sec_counter];

  PhaseSpace pp = Splitting(nBorn, COM);

  {
  // Test Phase space generation
  std::cout << "\033[1;37m--Test Phase space generation...\033[0m\n" << std::endl;
  clusterTree.print();
  pp.print();
  PhaseSpace pp_test = GenMomenta(pp, clusterTree);
  pp_test.print();


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

  {
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
  Vegas(3*(nBorn + nUnresolved) - 4, 1, *integrand_full2, &data2, 1, 0.01, 0.001, 0, 12, 100, sample_size, 100000, 10000, 1000, 2, "", &spin, &neval, &fail, integral, error, prob);
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
  }

  {
  // Test Infrared limits
  std::cout << "\033[1;37m--Test infrared limits...\033[0m\n" << std::endl;
  struct pair {
    TreeNode<Cluster>* node;
    Cluster data;
    bool operator < (const pair& pair2) const {
      return (this->data < pair2.data);
    }
  };
  std::cout << "\tLimiting variable eta:" << std::endl;
  for(int i = 0; i < nUnresolved; i++) {
    std::cout << std::setw(16) << "eta" << i << std::setw(9) << "eta" << i <<" (real)" << std::setw(17) << "Error" << " | ";
  }
  std::cout << std::endl;

  double scale = 1;
  double increment = 0.1;
  bool fail_eta = false;
  while (scale > 1.e-3) {
    scale *= increment;
    std::vector<std::vector<std::vector<double>>> xParFull;
    int level_int = 1;
    std::vector<TreeNode<Cluster>*> level = clusterTree.getLevel(level_int);
    std::vector<pair> pairs;
    for(TreeNode<Cluster>* node : level) {
      pair p;
      p.data = node->data;
      p.node = node;
      pairs.push_back(p);
    }
    std::sort(pairs.begin(), pairs.end());
    for(int i = 0; i < pairs.size(); i++) {
      level[i] = pairs[i].node;
    }
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
        if(node->children.size() > 1)
          xParFull.push_back(xPar);
      }
      level_int++;
      level = clusterTree.getLevel(level_int);
    }
    PhaseSpace ppFull = GenMomenta(pp, clusterTree, xParFull);
    level_int = 1;
    level = clusterTree.getLevel(level_int);
    while(level.size() > 0) {
      for(TreeNode<Cluster>* node : level){
        Momentum r = node->data.reference_momentum;
        for(Momentum u : node->data.unresolved_momenta) {
          double ratio = r*u/(r.components[0]*u.components[0]*2.)/std::pow(scale, level_int);
          std::cout << std::setw(17) << std::pow(scale, level_int) << std::setprecision(5)
                    << std::setw(17) << r*u/(r.components[0]*u.components[0]*2.) << std::setw(17) << std::abs(1 - ratio) <<  " | ";
          if(std::abs(ratio - 1) > 1.e-5) {
            fail_eta = true;
            std::cout << "fail_eta: " << std::abs(ratio - 1) << std::endl;
          }
        }
      }
      level_int++;
      level = clusterTree.getLevel(level_int);
    }
    std::cout << std::endl;
  }


  std::cout << "\n\tLimiting variable xi:" << std::endl;
  Momentum P = (-1.)*(pp.momenta[0] + pp.momenta[1]);
  for(int i = 0; i < nUnresolved; i++) {
    std::cout << std::setw(16) << "xi" << i << std::setw(9) << "xi" << i <<" (real)" << std::setw(17) << "Error" << " | ";
  }
  std::cout << std::endl;

  scale = 1;
  bool fail_xi = false;
  while (scale > 1.e-3) {
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
        if(node->children.size() > 1)
          xParFull.push_back(xPar);
      }
      level_int++;
      level = clusterTree.getLevel(level_int);
    }
    PhaseSpace ppFull = GenMomenta(pp, clusterTree, xParFull);

    level_int = 1;
    level = clusterTree.getLevel(level_int);

    while(level.size() > 0) {
      std::vector<pair> pairs;
      for(TreeNode<Cluster>* node : level) {
        pair p;
        p.data = node->data;
        p.node = node;
        pairs.push_back(p);
      }
      std::sort(pairs.begin(), pairs.end());
      for(int i = 0; i < pairs.size(); i++) {
        level[i] = pairs[i].node;
      }

      Momentum rTot, rWeighted, uTot;
      for(int j = 0; j < level.size(); j++){
        TreeNode<Cluster>* node = level[j];
        Momentum r = node->data.reference_momentum;
        rTot = rTot + r;
        for(int a = 0; a < node->data.unresolved; a++) {
          Momentum u = node->data.unresolved_momenta[a];
          uTot = uTot + u;
        }
        double xj = 2*r*(rWeighted - P)/(2.*P*(rWeighted - rTot - uTot) + (rTot + uTot)*(rTot + uTot) - rWeighted*rWeighted);
        for(int a = 0; a < node->data.unresolved; a++) {
          Momentum u = node->data.unresolved_momenta[a];
          Momentum urest;
          Momentum uhat = u/u.components[0];
          for(int k = 0; k < j; k++) for(int l = 0; l < level[k]->data.unresolved; l++) urest = urest + level[k]->data.unresolved_momenta[l];
          for(int l = 0; l < a; l++) urest = urest + level[j]->data.unresolved_momenta[l];

          double uMax = (2.*P*(rWeighted + r/xj - urest - rTot + r) + (urest + rTot - r)*(urest + rTot - r) - (rWeighted + r/xj)*(rWeighted + r/xj))/(2.*uhat*(P - rTot + r - urest));

          double ratio = u.components[0]/uMax/std::pow(scale, level_int);
          std::cout << std::setw(17) << std::pow(scale, level_int) << std::setprecision(5)
                    << std::setw(17) << u.components[0]/uMax << std::setw(17) << std::abs(1 - ratio) <<  " | ";
          if(std::abs(ratio - 1) > 1.e-5) {
            fail_xi = true;
            std::cout << "fail_xi: " << std::abs(ratio - 1) << std::endl;
          }
        }
        rWeighted = rWeighted + r/xj;
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
    std::cout << fail_xi << ", " << fail_eta << std::endl;
  }

  }
  }
  }
}