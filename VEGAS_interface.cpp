#include "VEGAS_interface.hpp"

int integrand_born(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata) {
  std::vector<double> xM, xcos, xphi;
  double *COM = static_cast<double*>(userdata);
  int n = (*ndim + 4)/3;
  for(int i = 0; i < *ndim; i++) {
    if(i < n - 2) {
      xM.push_back(x[i]);
    }
    else if ((i >= n - 2) and (i < 2*n - 3)) {
      xcos.push_back(x[i]);
    }
    else {
      xphi.push_back(x[i]);
    }
  }
  std::vector<std::vector<double>> xPar = {xM, xcos, xphi};
  PhaseSpace pp = Splitting(n, *COM, xPar);
  f[0] = pp.weight;
  return 0;
}

int integrand_full(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata) {
  std::vector<double> xM, xcos, xphi;
  UserData *data = static_cast<UserData*>(userdata);

  for(int i = 0; i < 3*(data->nBorn) - 4; i++) {
    if(i < data->nBorn - 2) {
      xM.push_back(x[i]);
    }
    else if ((i >= data->nBorn - 2) and (i < 2*data->nBorn - 3)) {
      xcos.push_back(x[i]);
    }
    else {
      xphi.push_back(x[i]);
    }
  }
  std::vector<std::vector<double>> xPar = {xM, xcos, xphi};
  PhaseSpace pp = Splitting(data->nBorn, data->COM, xPar);

  int counter = 3*(data->nBorn) - 4;
  std::vector<std::vector<std::vector<double>>> xPar_rest;
  for(int i = 0; i < data->cluster.size(); i++) {
    std::vector<std::vector<double>> xPar_cluster;
    for(int j = 0; j < data->cluster[i].unresolved; j++) {
      std::vector<double> xPar_unresolved(3);
      xPar_unresolved[0] = x[counter];
      xPar_unresolved[1] = x[counter + 1];
      xPar_unresolved[2] = x[counter + 2];
      counter += 3;
      xPar_cluster.push_back(xPar_unresolved);
    }
    xPar_rest.push_back(xPar_cluster);
  }
  PhaseSpace pp_full = GenMomenta(pp, data->cluster, xPar_rest);
  f[0] = pp_full.weight;
  return 0;
}

int integrand_full2(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata) {
  std::vector<double> xM, xcos, xphi;
  UserData2 *data = static_cast<UserData2*>(userdata);

  for(int i = 0; i < 3*(data->nBorn) - 4; i++) {
    if(i < data->nBorn - 2) {
      xM.push_back(x[i]);
    }
    else if ((i >= data->nBorn - 2) and (i < 2*data->nBorn - 3)) {
      xcos.push_back(x[i]);
    }
    else {
      xphi.push_back(x[i]);
    }
  }
  std::vector<std::vector<double>> xPar = {xM, xcos, xphi};
  PhaseSpace pp = Splitting(data->nBorn, data->COM, xPar);

  int counter = 3*(data->nBorn) - 4;
  std::vector<std::vector<std::vector<double>>> xPar_rest;
  std::vector<TreeNode<Cluster>*> nodes = data->clusterTree->getNodes();
  std::vector<Cluster> cluster;
  for(auto& node : nodes) cluster.push_back(node->data);
  for(int j = 0; j < cluster.size(); j++) {
    std::vector<std::vector<double>> xPar_cluster;
    for(int a = 0; a < cluster[j].unresolved; a++) {
      std::vector<double> xPar_unresolved(3);
      xPar_unresolved[0] = x[counter];
      xPar_unresolved[1] = x[counter + 1];
      xPar_unresolved[2] = x[counter + 2];
      counter += 3;
      xPar_cluster.push_back(xPar_unresolved);
    }
    xPar_rest.push_back(xPar_cluster);
  }
  PhaseSpace pp_full = GenMomenta2(pp, *(data->clusterTree), xPar_rest);
  f[0] = pp_full.weight;
  return 0;
}