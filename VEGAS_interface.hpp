#ifndef VEGAS_INTERFACE_HPP
#define VEGAS_INTERFACE_HPP

#include <iostream>
#include "src/PhaseSpace.hpp"
#include "src/Tree.hpp"
#include "/home/tom/Documents/software/software/Cuba-4.2.2/cuba.h"

using namespace PSF;

struct UserData {
  double COM;
  std::vector<Cluster> cluster;
  int nBorn;
  UserData(double COM, int nBorn, std::vector<Cluster> cluster):COM(COM), nBorn(nBorn), cluster(cluster){}
  UserData(){};
};

struct UserData2 {
  double COM;
  Tree<Cluster> *clusterTree;
  int nBorn;
  int nUnresolved;
  UserData2(double COM, int nBorn, int nUnresolved, Tree<Cluster> *clusterTree):COM(COM), nBorn(nBorn), nUnresolved(nUnresolved), clusterTree(clusterTree){}
  UserData2(){};
};

int integrand_born(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata);

int integrand_full(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata);

int integrand_full2(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata);

#endif
