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

  // Generate custom phase space point with Tree
  std::cout << "\n################################################################\n" << std::endl;
  std::cout << "Tree (Manually)\n" << std::endl;
  std::cout << "################################################################\n" << std::endl;
  // This phase-space generation algorithm takes a Born phase-space point (generated e.g. with RAMBO
  // or the above Splitting algorithm) and adds a custom number of additional partons. The phase-space
  // is uses the (xi, eta) variables from the Sector improved residue subtraction scheme but applies
  // them recursively to Clusters of partons.

  // Example 1: 4 outgoing momenta in Born configuration + 1 additional radiation
  int nBorn = 4;
  int nUnresolved = 1;
  std::cout << "Example: 2 ->  " << nBorn << " + " << nUnresolved << std::endl;
  std::cout << "Born configuration phase-space point:" << std::endl;
  PhaseSpace pp = RAMBO(nBorn, COM);
  pp.print();

  // The Root node does not yet hold real information!
  std::unique_ptr<TreeNode<Cluster>> root =  std::make_unique<TreeNode<Cluster>>(Cluster());
  root->data.isReference = true;
  std::unique_ptr<TreeNode<Cluster>> node0 =  std::make_unique<TreeNode<Cluster>>(Cluster(3, 1)); // First index is the index of the reference parton. Second index is the number of unresolved Clusters
  node0->data.isReference = true;
  std::unique_ptr<TreeNode<Cluster>> node1 = std::make_unique<TreeNode<Cluster>>(Cluster(3, 1));
  node1->data.isReference = true;
  std::unique_ptr<TreeNode<Cluster>> node2 = std::make_unique<TreeNode<Cluster>>(Cluster(4, 1));

  // Create tree from nodes
  Tree<Cluster> clusterTreeCustom;
  clusterTreeCustom.setRoot(std::move(root));
  clusterTreeCustom.addChild(clusterTreeCustom.getRoot(), std::move(node0));
  clusterTreeCustom.addChild(clusterTreeCustom.getRoot()->children[0].get(), std::move(node1));
  clusterTreeCustom.addChild(clusterTreeCustom.getRoot()->children[0].get(), std::move(node2));
  std::cout << "\nTree of Clusters used for the generation of the Full phase-space point" << std::endl;
  clusterTreeCustom.print(); // Print resulting tree

  // Define integration parameters (eta, xi, phi)
  // First index labels the level
  // Second index labes the unresolved Cluster in the level
  // Third index labels eta, xi , phi
  std::vector<std::vector<std::vector<double>>> xParFull;
  int level_int = 1;
  std::vector<PSF::TreeNode<PSF::Cluster>*> level = clusterTreeCustom.getLevel(level_int);
  int unresolved_counter = 0;
  while(level.size() > 0) {
    int unresolved_level = 0;
    for(PSF::TreeNode<PSF::Cluster>* node : level){
      unresolved_level += node->data.unresolved;

      std::vector<std::vector<double>> xPar;
      for(int c = 0; c < node->data.unresolved; c++) {
        double xi, eta, phi;
        // Specify xi, eta, phi values. All variables are in [0, 1]
        xi = 0.5; // Energy variable
        eta = 0.01; // Azimuthal angle variable
        phi = 0.4; // Polar angle variable
        unresolved_counter++;
        std::vector<double> xPar_c = {eta, xi, phi};
        xPar.push_back(xPar_c);
      }
      if(node->children.size() > 1)
        xParFull.push_back(xPar);
    }
    level_int++;
    // update level
    level = clusterTreeCustom.getLevel(level_int);
  }

  PSF::PhaseSpace ppFull = PSF::GenMomenta(pp, clusterTreeCustom, xParFull);
  std::cout << "\nFull phase-space point: " << std::endl;
  ppFull.print();

  // Generate custom phase space point with Tree
  std::cout << "\n################################################################\n" << std::endl;
  std::cout << "Tree (Automatic)\n" << std::endl;
  std::cout << "################################################################\n" << std::endl;
  // This phase-space generation algorithm takes a Born phase-space point (generated e.g. with RAMBO
  // or the above Splitting algorithm) and adds a custom number of additional partons. The phase-space
  // is uses the (xi, eta) variables from the Sector improved residue subtraction scheme but applies
  // them recursively to Clusters of partons.
  std::cout << "\nThe trees can also be generated automatically.\n" << std::endl;
  nUnresolved = 2;
  std::cout << "Example: 2 ->  " << nBorn - 2 << " + " << nUnresolved << std::endl;
  // Define ClusterTree
  const std::vector<bool> flavor = {0,0,1,1,1,1}; // this vector tells the tree generation which of the outgoing partons are allowed to split. 0 = Not allowed, 1 = allowed
  std::vector<PSF::Tree<PSF::Cluster>> trees = PSF::GenTrees(nUnresolved);
  // print all possible trees

  for(int tree_counter = 0; tree_counter < trees.size(); tree_counter++) {
    PSF::Tree<PSF::Cluster> tree = trees[tree_counter];
    std::cout << "tree_counter = " << tree_counter << std:: endl;
    tree.print();

    // The tree object is bare in the sense that the momenta are not yet assigned to the nodes. This is done in GenSectors.
    std::vector<PSF::Tree<PSF::Cluster>> sectors = PSF::GenSectors(flavor, tree, nBorn + 2);

    int sec_counter = 0;
    for(int sec_counter = 0; sec_counter < sectors.size(); sec_counter++) {
      PSF::Tree<PSF::Cluster> clusterTree = sectors[sec_counter];

      xParFull.clear();
      level_int = 1;
      level = clusterTree.getLevel(level_int);
      unresolved_counter = 0;
      while(level.size() > 0) {
        int unresolved_level = 0;
        for(PSF::TreeNode<PSF::Cluster>* node : level){
          unresolved_level += node->data.unresolved;

          std::vector<std::vector<double>> xPar;
          for(int c = 0; c < node->data.unresolved; c++) {
            double xi, eta, phi;
            // Specify xi, eta, phi values. All variables are in [0, 1]
            xi = 0.5; // Energy variable
            eta = 0.001; // Azimuthal angle variable
            phi = 0; // Polar angle variable
            unresolved_counter++;
            std::vector<double> xPar_c = {eta, xi, phi};
            xPar.push_back(xPar_c);
          }
          if(node->children.size() > 1)
            xParFull.push_back(xPar);
        }
        level_int++;
        // update level
        level = clusterTree.getLevel(level_int);
      }
      ppFull = PSF::GenMomenta(pp, clusterTree, xParFull);
      std::cout << "\nFull phase-space point: " << std::endl;
      ppFull.print();
      std::cout << "\nCorresponding weight: " << ppFull.weight << std::endl;
    }
  }
}