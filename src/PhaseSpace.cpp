#include "PhaseSpace.hpp"

void LorentzMatrix::print() const {
  int precision = 6;
  std::cout << std::setprecision(precision);
  for(int i = 0; i < 4; i++) {
  if(i == 0) std::cout << "/";
  else if (i == 4 - 1) std::cout << "\\";
  else std::cout << "|" ;
  for(int j = 0; j < 4; j++) {
    std::cout << std::setw(precision+5) << components[i][j] ;
    if(j != 4 - 1) std::cout << " ";
  }
  if(i == 0) std::cout << "\\" << std::endl;
  else if (i == 4 - 1) std::cout << "/" << std::endl;
  else std::cout << "|" << std::endl;
  }
}

void LorentzMatrix::print_list() const {
  std::cout << "{";
  for(int i = 0; i < 4; i++) {
    std::cout << "{";
    for(int j = 0; j < 4; j++) {
      std::cout << components[i][j];
      if(j != 4 - 1) std::cout << ", ";
      else std::cout << "}";
    }
    if(i != 4 - 1) std::cout << ", ";
  }
  std::cout << "}" << std::endl;
}

void PhaseSpace::check_momentum_conservation(double acc) {
  Momentum check;
  for(int i = 0; i < size(); i++) {
    check = check + momenta[i];
  }
  for(int i = 0; i < 4; i++) {
    if(check.components[i] > acc) {
      std::cout << "WARNING: Momentum not conserved p = ";
      check.print();
      return;
    }
  }
  std::cout << "Momentum is conserved p = ";
  check.print();
}

void PhaseSpace::check_onshellness(double acc) {
  for(int i = 0; i < size(); i++) {
    std::cout << "p[" << i << "]^2 = " << momenta[i].inv_mass2() << std::endl;
  }
}

void PhaseSpace::print() const {
  for(int i = 0; i < size(); i++) {
    std::cout << "p[" << i << "] = ";
    momenta[i].print();
  }
}

Momentum operator+(Momentum p1, Momentum p2) {
  return Momentum(std::vector<double>({p1.components[0] + p2.components[0], p1.components[1] + p2.components[1], p1.components[2] + p2.components[2], p1.components[3] + p2.components[3]}));
}

Momentum operator-(Momentum p1, Momentum p2) {
  return Momentum(std::vector<double>({p1.components[0] - p2.components[0], p1.components[1] - p2.components[1], p1.components[2] - p2.components[2], p1.components[3] - p2.components[3]}));
}

double operator*(Momentum p1, Momentum p2) {
  return p1.components[0]*p2.components[0] - p1.components[1]*p2.components[1] - p1.components[2]*p2.components[2] - p1.components[3]*p2.components[3];
}

Momentum operator*(double a, Momentum p) {
  return Momentum(std::vector<double>({a*p.components[0], a*p.components[1], a*p.components[2], a*p.components[3]}));
}

Momentum operator*(Momentum p, double a) {
  return a*p;
}

Momentum operator-(Momentum p) {
  return (-1.)*p;
}

Momentum operator/(Momentum p, double a) {
  return (1./a)*p;
}

Momentum operator*(LorentzMatrix lam, Momentum p) {
  Momentum p_out;
  for(int i = 0; i < 4; i++) for(int dummy = 0; dummy < 4; dummy++) {
    p_out.components[i] += lam.components[i][dummy]*p.components[dummy];
  }
  return p_out;
}

LorentzMatrix operator*(LorentzMatrix lam1, LorentzMatrix lam2) {
  LorentzMatrix lam_out;
  for(int i = 0; i < 4; i++) for(int j = 0; j < 4; j++) for(int dummy = 0; dummy < 4; dummy++) {
    lam_out.components[i][j] += lam1.components[i][dummy]*lam2.components[dummy][j];
  }
  return lam_out;
}

double rnd(double lower, double upper) {
  double f = (double)rand() /RAND_MAX;
  return lower + f *(upper - lower);
}

double RAMBO_measure(int nMomenta, double COM) {
  double volume = std::pow(M_PI/2., nMomenta - 1)/double(factorial(nMomenta - 1))/double(factorial(nMomenta - 2))*std::pow(COM, 2*nMomenta - 4); //Byckling and Kajantie (VI.2.17)
  double normalization = std::pow(2.*M_PI, 4)/std::pow(2.*M_PI, 3*nMomenta);
  return volume*normalization;
}

PhaseSpace RAMBO(int nMomenta, double COM) {
  // generate n massless momenta uniformly distributed with energy following E*Exp(-E)
  std::vector<Momentum> massless_momenta;

  for(int i = 0; i < nMomenta; i++) {
    double rho1 = rnd(0., 1.);
    double rho2 = rnd(0., 1.);
    double rho3 = rnd(0., 1.);
    double rho4 = rnd(0., 1.);

    double cos = 2*rho1 - 1.;
    double sin = std::sqrt(1. - cos*cos);
    double phi = 2*M_PI*rho2;
    double q0 = -std::log(rho3*rho4);
    double qx = q0*sin*std::cos(phi);
    double qy = q0*sin*std::sin(phi);
    double qz = q0*cos;
    massless_momenta.push_back(Momentum(std::vector<double>({q0, qx, qy, qz})));
  }
  Momentum Q;
  for(int i = 0; i < nMomenta; i++) {
    Q = Q + massless_momenta[i];
  }
  double M = Q.inv_mass();
  double x = COM/M;
  Momentum b = -1.*Q/M; b.components[0] = 0;
  double gamma = Q.components[0]/M;
  double a = 1./(1. + gamma);
  std::vector<Momentum> output;
  output.push_back(Momentum(std::vector<double>({-COM/2., 0, 0, -COM/2.})));
  output.push_back(Momentum(std::vector<double>({-COM/2., 0, 0, COM/2.})));

  for(int i = 0; i < nMomenta; i++) {
    std::vector<double> pi = {0,0,0,0};
    pi[0] = x*(gamma*massless_momenta[i].components[0] - b*massless_momenta[i]);
    for(int j = 1; j < 4; j++) {
      pi[j] = x*(massless_momenta[i].components[j] + massless_momenta[i].components[0]*b.components[j] - a*(b*massless_momenta[i])*b.components[j]);
    }
    output.push_back(Momentum(pi));
  }
  PhaseSpace pp(output);
  pp.weight = RAMBO_measure(nMomenta, COM);
  return pp;
}

// Lam becomes rotation matrix
LorentzMatrix rotation(double x, double y, double z) {
  LorentzMatrix Lam;
  Lam.components[0][0] = 1; Lam.components[0][1] = 0; Lam.components[0][2] = 0; Lam.components[0][3] = 0;
  Lam.components[1][0] = 0; Lam.components[1][1] = std::cos(z)*std::cos(y);  Lam.components[1][2] = std::cos(z)*std::sin(y)*std::sin(x)-std::sin(z)*std::cos(x);  Lam.components[1][3] = std::cos(z)*std::sin(y)*std::cos(x)+std::sin(z)*std::sin(x);
  Lam.components[2][0] = 0; Lam.components[2][1] = std::sin(z)*std::cos(y); Lam.components[2][2] = std::sin(z)*std::sin(y)*std::sin(x)+std::cos(z)*std::cos(x); Lam.components[2][3] = std::sin(z)*std::sin(y)*std::cos(x)-std::cos(z)*std::sin(x);
  Lam.components[3][0] = 0; Lam.components[3][1] = -std::sin(y); Lam.components[3][2] = std::cos(y)*std::sin(x); Lam.components[3][3] = std::cos(y)*std::cos(x);
  return Lam;
}

LorentzMatrix find_LT(Momentum v1, Momentum v2) {
  // Returns Lambda to be the Lorentz transformation matrix between v1 = Lambda * v2
  //First check if vectors are compatible
  double m1 = v1.inv_mass();
  double m2 = v2.inv_mass();
  if ( std::abs(m1 - m2) > 3.e-5) {
    std::cout << "Warning: Momenta are not compatible" << std::endl;
    std::cout << "m1 = \t" << m1 << "\t m2 = \t" << m2 << std::endl;
    std::cout.precision(11);
    std::cout << "v1 = (" << v1.components[0] << ", " << v1.components[1] << ", " << v1.components[2] << ", " << v1.components[3] << ")" << std::endl;
    std::cout << "v2 = (" << v2.components[0] << ", " << v2.components[1] << ", " << v2.components[2] << ", " << v2.components[3] << ")" << std::endl;
  }
  //Next we determine the rotational part of the transformation
  //We determine the rotational axis by doing a cross product
  LorentzMatrix R;
  Momentum v1_hat, v2_hat;
  double v1_abs = std::sqrt(v1.components[1] * v1.components[1] + v1.components[2] * v1.components[2] + v1.components[3] * v1.components[3]);
  double v2_abs = std::sqrt(v2.components[1] * v2.components[1] + v2.components[2] * v2.components[2] + v2.components[3] * v2.components[3]);
  if (v1_abs > 1.E-12 and v2_abs > 1.E-12) {
    double phi_2 = std::atan2(v2.components[2],v2.components[1]);
    double theta_2 = std::acos(v2.components[3]/v2_abs);
    double theta_1 = std::acos(v1.components[3]/v1_abs);
    double phi_1 = std::atan2(v1.components[2],v1.components[1]);

    LorentzMatrix Lam_help, Lam_help2;
    R = rotation(0, 0, -phi_2);
    Lam_help = rotation(0, theta_1-theta_2, 0);
    Lam_help2 = Lam_help*R;
    Lam_help = rotation(0, 0, phi_1);
    R = Lam_help*Lam_help2;

    v1_hat = v1/ v1_abs; //Normalize v1_hat
    v2_hat = v2/ v2_abs; //Normalize v2_hat
  }
  else {
    //std::cout << "Fail safe solution (v1, v2 < 1.E-10); v1_abs = \t" << v1_abs << "\t v2_abs = \t" << v2_abs << std::endl;
    if (v1_abs < 1.E-10 and v2_abs > 1.E-10) {
      v2_hat = v2/ v2_abs; //Normalize v2_hat
    }
    else if (v1_abs > 1.E-10 and v2_abs < 1.E-10) {
      v1_hat = v1/ v1_abs; //Normalize v1_hat
    }
    for(int i = 0; i < 4; i++) for(int j = 0; j < 4; j++) {
      if(i == j) R.components[i][j] = 1.;
      else R.components[i][j] = 0.;
    }
  }

  //Next we work out the boost part
  double m = (m1 + m2)/2; //average the mass which should be the same anyway
  //double beta = (- v2.components[0] * v2_abs + v1.components[0] * v1_abs ) / (m * m - v1.components[0] * v1.components[0]  - v2.components[0] * v2.components[0] );    //Lorentz beta
  double beta = -(v1.components[0] * v1_abs - v2.components[0] * v2_abs )/(v1_abs * v1_abs + v2.components[0] * v2.components[0]);
  double gamma = 1/std::sqrt(1 - beta * beta);   //Lorentz gamma
  double limit = (gamma - 1.)/beta/beta;
  if (std::abs(beta) < 1.e-7) {
      limit = 1./2.;
  }
  // determine boost vector which alligns with v1_hat since this is the reference vector
  Momentum beta_vec;
  if (v1_abs < 1.E-12 and v2_abs > 1.E-12) {
    beta_vec = v2_hat*beta;
  }
  else if (v1_abs > 1.E-12 and v2_abs < 1.E-12) {
    beta_vec = v1_hat*beta;
  }
  else if (v1_abs < 1.E-12 and v2_abs < 1.E-12) {
    beta_vec = v2_hat*beta;
    std::cout << "Lorentz transformation is identity" << std::endl;
  }
  else {
    beta_vec=v1_hat*beta;
  }

  LorentzMatrix boost(std::vector<std::vector<double>>({{gamma, -gamma * beta_vec.components[1], -gamma * beta_vec.components[2], -gamma * beta_vec.components[3]},
    {-gamma * beta_vec.components[1], 1 + limit * beta_vec.components[1] * beta_vec.components[1], limit * beta_vec.components[1] * beta_vec.components[2], limit * beta_vec.components[1] * beta_vec.components[3]},
    {-gamma * beta_vec.components[2], limit * beta_vec.components[2] * beta_vec.components[1], 1 + limit * beta_vec.components[2] * beta_vec.components[2], limit * beta_vec.components[2] * beta_vec.components[3]},
    {-gamma * beta_vec.components[3], limit * beta_vec.components[3] * beta_vec.components[1], limit * beta_vec.components[3] * beta_vec.components[2], 1 + limit * beta_vec.components[3] * beta_vec.components[3]}}));

  // Combine boost and rotation
  LorentzMatrix Lambda = boost*R;
  return Lambda;
}

PhaseSpace Splitting(int nMomenta, double COM, std::vector<std::vector<double>> x) {
  // Generate Phase space point with nMomenta momenta
  // x must be of length {{n - 2}, {n - 1}, {n - 1}}
  PhaseSpace output;

  std::vector<double> M(nMomenta - 1);
  std::vector<double> cos(nMomenta - 1);
  std::vector<double> phi(nMomenta - 1);

  M[0] = COM;
  for(int i = 1; i < nMomenta - 1; i++) {
    M[i] = x[0][i-1]*M[i-1];
  }
  for(int i = 0; i < x[1].size(); i++) cos[i] = 1. - 2.*x[1][i];
  for(int i = 0; i < x[2].size(); i++) phi[i] = 2.*M_PI*x[2][i];

  std::reverse(M.begin(), M.end());
  std::reverse(cos.begin(), cos.end());
  std::reverse(phi.begin(), phi.end());

  std::vector<Momentum> momenta;
  double jacobian = 1./2./COM;
  for(int i = 0; i < nMomenta - 1; i++) {
    double sin = std::sqrt(1 - std::pow(cos[i],2));
    double p0;
    if(i == 0) p0 = M[0]/2.;
    else {
      p0 = (std::pow(M[i], 2) - std::pow(M[i - 1], 2))/2./M[i];
      jacobian *= M[i];
    }
    jacobian *= p0/2.;
    jacobian *= 2.*2.*M_PI;
    double px = p0*sin*std::cos(phi[i]);
    double py = p0*sin*std::sin(phi[i]);
    double pz = p0*cos[i];
    Momentum p(std::vector<double>({p0,px,py,pz}));

    if(i != 0) {
      Momentum pCOM(std::vector<double>({M[i], 0, 0, 0}));
      Momentum pRest = pCOM - p;
      Momentum pCOM_old = Momentum(std::vector<double>({M[i - 1], 0, 0, 0}));
      // Find Transformation yielding pRest from (M_{i-1},0)
      LorentzMatrix lam = find_LT(pRest, pCOM_old);
      for(int j = 0; j < momenta.size(); j++) {
        momenta[j] = lam*momenta[j];
      }
      momenta.push_back(p);
    }
    else {
      momenta.push_back(p);
      momenta.push_back(Momentum(std::vector<double>({p0, -px, -py, -pz})));
    }
  }
  momenta.insert(momenta.begin(), Momentum(std::vector<double>({-COM/2., 0, 0, COM/2.})));
  momenta.insert(momenta.begin(), Momentum(std::vector<double>({-COM/2., 0, 0, -COM/2.})));

  output.momenta = momenta;
  double normalization = std::pow(2.*M_PI, 4)/std::pow(2.*M_PI, 3*nMomenta);
  output.weight = normalization*jacobian;
  return output;
}

PhaseSpace Splitting(int nMomenta, double COM) {
  std::vector<double> xM, xcos, xphi;

  for(int i = 0; i < 3*nMomenta - 4; i++) {
    if(i < nMomenta - 2) {
      xM.push_back(rnd(0, 1));
    }
    else if ((i >= nMomenta - 2) and (i < 2*nMomenta - 3)) {
      xcos.push_back(rnd(0, 1));
    }
    else {
      xphi.push_back(rnd(0, 1));
    }
  }
  std::vector<std::vector<double>> xPar = {xM, xcos, xphi};
  return Splitting(nMomenta, COM, xPar);
}

PESCPhaseSpace GenMomenta(const PhaseSpace pp, const std::vector<Cluster>& cluster, std::vector<std::vector<std::vector<double>>> xPar) {
  PESCPhaseSpace pp_new;
  pp_new.cluster = cluster;
  for(int j = 0; j < pp_new.cluster.size(); j++) {
    pp_new.cluster[j].unresolved_momenta = std::vector<Momentum>(0);
  }
  std::vector<double> x(pp_new.cluster.size(), 0.);
  int dim = 4;
  double jacobian = 1.;
  int nUnresolved = 0;
  Momentum P = (-1.)*(pp.momenta[0] + pp.momenta[1]);
  Momentum rTot, rWeighted, uTot;
  for(int j = 0; j < pp_new.cluster.size(); j++) {
    Momentum r = pp.momenta[pp_new.cluster[j].reference];
    // Find rotation that yields r from (1, 0, 0, 1) and apply to uhat
    LorentzMatrix rotation = find_LT(r/r.components[0], Momentum(std::vector<double>({1., 0., 0., 1.})));
    rTot = rTot + r;
    // Generate the unresolved momenta in the cluster
    for(int a = 0; a < pp_new.cluster[j].unresolved; a++) {
      double eta = xPar[j][a][0];
      double xi = xPar[j][a][1];
      double phi = xPar[j][a][2]*2.*M_PI;
      double cos = 1. - 2.*eta;
      double sin = std::sqrt(1. - cos*cos);
      Momentum uhat(std::vector<double>({1., sin*std::cos(phi), sin*std::sin(phi), cos}));
      uhat = rotation*uhat;

      Momentum urest;
      for(int k = 0; k < j; k++) for(int l = 0; l < pp_new.cluster[k].unresolved; l++) urest = urest + pp_new.cluster[k].unresolved_momenta[l];
      for(int l = 0; l < a; l++) urest = urest + pp_new.cluster[j].unresolved_momenta[l];

      double uMax = (2.*P*rTot - rTot*rTot + urest*urest - 2.*P*urest - 2.*P*rWeighted + rWeighted*rWeighted + 2.*rWeighted*urest)/
          (2.*uhat*(P - rWeighted - urest));

      double u0 = uMax*xi;
      Momentum u = u0*uhat;
      jacobian *= u0*u0*2.*2.*M_PI*uMax /2./std::pow(2.*M_PI, 3)/u0;
      pp_new.cluster[j].unresolved_momenta.push_back(u);
      nUnresolved++;
      uTot = uTot + u;
    }

    double xj = (-2.*P*rTot + rTot*rTot - uTot*uTot + 2.*P*uTot + 2.*P*rWeighted - rWeighted*rWeighted - 2.*rWeighted*uTot)/
        (-2.*P*r + 2.*r*rWeighted + 2.*r*uTot);
    x[j] = xj;

    // determine weight factor
    double c1 = -2.*P*r + 2.*r*rTot;
    double c2 = -2.*P*r + 2.*r*rWeighted + 2.*r*uTot;
    double c3 = -2.*P*rWeighted + rWeighted*rWeighted - 2.*P*uTot + uTot*uTot + 2.*rWeighted*uTot + 2.*P*(rTot - r) - (rTot - r)*(rTot - r);

    jacobian *= pow(xj, dim - 3)* (-1.)/(-c1*c2/pow(c2*xj + c3, 2));
    if(xj < 0) {
      jacobian = 0.;
      std::cout << "xj = " << xj << std::endl;
    }
    rWeighted = rWeighted + r*x[j];
    pp_new.cluster[j].reference_momentum = x[j]*r;
  }

  // determine q and qTilde
  Momentum qFull = P, q = P;
  for(int j = 0; j < pp_new.cluster.size(); j++) {
    qFull = qFull - pp_new.cluster[j].reference_momentum;
    q = q - pp.momenta[pp_new.cluster[j].reference];
    for(int k = 0; k < pp_new.cluster[j].unresolved; k++) {
      qFull = qFull - pp_new.cluster[j].unresolved_momenta[k];
    }
  }

  // determine Lorentz transformation yielding qFull from q
  LorentzMatrix Lambda = find_LT(qFull, q);
  pp_new.momenta = pp.momenta;
  pp_new.momenta[0] = pp.momenta[0];
  pp_new.momenta[1] = pp.momenta[1];

  // Apply Lorentz transformation to all final states momenta that are not reference vectors
  std::vector<int> qi_indices;
  for(int i = 2; i < pp.size(); i++) {
    bool reference = false;
    int cluster_index;
    for(int j = 0; j < pp_new.cluster.size(); j++){
      if(pp_new.cluster[j].reference == i) {
        reference = true;
        cluster_index = j;
      }
    }
    if(!reference) {
      Momentum pi_new = Lambda*pp.momenta[i];
      pp_new.momenta[i] = pi_new;
    }
    else {
      pp_new.momenta[i] = pp_new.cluster[cluster_index].reference_momentum;
    }
  }
  for(int i = 0; i < pp_new.cluster.size(); i++) {
    for(Momentum& u : pp_new.cluster[i].unresolved_momenta) {
      pp_new.momenta.push_back(u);
    }
  }
  pp_new.weight = pp.weight*jacobian;
  return pp_new;
}

PESCPhaseSpace GenMomenta(const PhaseSpace pp, const std::vector<Cluster>& cluster) {
  std::vector<std::vector<std::vector<double>>> xPar;
  for(int j = 0; j < cluster.size(); j++) {
    std::vector<std::vector<double>> xPar_Cluster;
    for(int a = 0; a < cluster[j].unresolved; a++) {
      std::vector<double> xPar_unresolved(3);
      xPar_unresolved[0] = rnd(0, 1);
      xPar_unresolved[1] = rnd(0, 1);
      xPar_unresolved[2] = rnd(0, 1);
      xPar_Cluster.push_back(xPar_unresolved);
    }
    xPar.push_back(xPar_Cluster);
  }
  return GenMomenta(pp, cluster, xPar);
}

PhaseSpace GenMomenta2(const PhaseSpace pp, const Tree<Cluster>& clusterTree, std::vector<std::vector<std::vector<double>>> xPar) {
  TreeNode<Cluster>* root = clusterTree.getRoot();
  PhaseSpace pp_new = pp;
  TreeNode<Cluster>* current = root;
  int level_int = 1;
  int clusterCounter = 0;
  std::vector<TreeNode<Cluster>*> level = clusterTree.getLevel(level_int);
  int unresolved_counter = pp.size();
  int last_unresolved_counter = unresolved_counter;
  while(level.size() > 0) {
    std::vector<Cluster> currCluster;
    for(TreeNode<Cluster>* cluster : level) {
      if(cluster->children.size() > 1)
        currCluster.push_back(cluster->data);
    }
    std::vector<std::vector<std::vector<double>>> currxPar;
    for(int j = 0; j < currCluster.size(); j++) {
      if(currCluster[j].unresolved > 0) {
        currxPar.push_back(xPar[clusterCounter]);
        clusterCounter++;
      }
    }

    if(currCluster.size() != 0) {
      pp_new = GenMomenta(pp_new, currCluster, currxPar);
    }
    for(TreeNode<Cluster>* cluster : level) {
      cluster->data.reference_momentum = pp_new.momenta[cluster->data.reference];
      cluster->data.unresolved_momenta = std::vector<Momentum>(cluster->data.unresolved);
      while(unresolved_counter < last_unresolved_counter + cluster->data.unresolved) {
        cluster->data.unresolved_momenta[unresolved_counter - last_unresolved_counter] = pp_new.momenta[unresolved_counter];
        unresolved_counter++;
      }
      last_unresolved_counter = unresolved_counter;
    }

    level_int++;
    level = clusterTree.getLevel(level_int);
  }
  return pp_new;
}

PhaseSpace GenMomenta2(const PhaseSpace pp, const Tree<Cluster>& clusterTree) {
  std::vector<std::vector<std::vector<double>>> xPar;
  std::vector<TreeNode<Cluster>*> nodes = clusterTree.getNodes();
  std::vector<Cluster> cluster;
  for(auto& node : nodes) {
    if(node->children.size() > 0)
      cluster.push_back(node->data);
  }
  for(int j = 0; j < cluster.size(); j++) {
    std::vector<std::vector<double>> xPar_Cluster;
    for(int a = 0; a < cluster[j].unresolved; a++) {
      std::vector<double> xPar_unresolved(3);
      xPar_unresolved[0] = rnd(0, 1);
      xPar_unresolved[1] = rnd(0, 1);
      xPar_unresolved[2] = rnd(0, 1);
      xPar_Cluster.push_back(xPar_unresolved);
    }
    if(cluster[j].unresolved > 0)
      xPar.push_back(xPar_Cluster);
  }

  return GenMomenta2(pp, clusterTree, xPar);
}

std::vector<Tree<Cluster>> GenTrees(const Tree<Cluster>& tree, int nUnresolved) {
  // Based on the input "tree", this function generates all possible next levels with <=nUnresolved unresolved nodes.
  std::vector<Tree<Cluster>> output;
  if(nUnresolved <= 0) return output;
  std::vector<TreeNode<Cluster>*> level = tree.getLevel(tree.depth() - 1);
  // Find out wich of the last nodes are reference
  std::vector<bool> isReference;
  int nUnresolved_sofar = 0;
  for(TreeNode<Cluster>* node : level) {
    if(node->data.isReference)
      isReference.push_back(true);
    else {
      isReference.push_back(false);
    }
  }
  for(TreeNode<Cluster>* node : tree.getNodes()) {
    if(node == tree.getRoot()) continue;
    if(!node->data.isReference)
      nUnresolved_sofar++;
  }

  // Every unresolved can split, ergo the max number of new nodes is 2*nUresolved
  for(int i = 1; i <= 2*nUnresolved; i++) {
    // A patition of i corresponds to the way we distribute the i unresolved nodes among the existing nodes
    std::vector<std::vector<int>> partitions = getPartitions(i);
    for(std::vector<int> partition : partitions) {
      if(partition.size() > level.size()) continue; // Skip partitions with more nodes than existing in the level
      // Fill remaining spots with zero
      for(int dummy = partition.size(); dummy < level.size(); dummy++) partition.push_back(0);
      // The partition can appear in all possible permutations
      std::vector<std::vector<int>> permutations = getPermutations(partition);
      for(std::vector<int>& permutation : permutations) {
        bool valid_splitting = true;
        Tree<Cluster> treeCopy(tree);
        for(int node_counter = 0; node_counter < level.size(); node_counter++) {
          TreeNode<Cluster>* node = treeCopy.getLevel(tree.depth() - 1)[node_counter];
          bool generate_unresolved = false;
          // If a reference node splits, we need a new reference node in the next level
          if(isReference[node_counter] and (permutation[node_counter] > 0)) {
            TreeNode<Cluster>* child_ref = new TreeNode<Cluster>(true);
            treeCopy.addChild(node, child_ref);
            generate_unresolved = true;
          }
          // Remove trees where a reference has never any splitting
          if(isReference[node_counter] and (node->parent == treeCopy.getRoot()) and (permutation[node_counter] == 0)) {
            valid_splitting = false;
            break;
          }
          // Remove trees where an unresolved does not split
          if(!isReference[node_counter] and (permutation[node_counter] == 1)) {
            valid_splitting = false;
            break;
          }
          // Generate unresolved nodes in the next level
          for(int n = 0; n < permutation[node_counter]; n++) {
            if(!generate_unresolved) {
              TreeNode<Cluster>* child_ref = new TreeNode<Cluster>(true);
              treeCopy.addChild(node, child_ref);
              generate_unresolved = true;
            }
            else {
              TreeNode<Cluster>* child = new TreeNode<Cluster>(false);
              treeCopy.addChild(node, child);
            }
          }
        }
        int nUnresolved_new = 0;
        for(TreeNode<Cluster>* node : treeCopy.getNodes()) {
          if(node != treeCopy.getRoot()) {
            if(!node->data.isReference) nUnresolved_new++;
          }
        }
        // Remove trees with two many unresolved nodes
        if(nUnresolved_new - nUnresolved_sofar > nUnresolved) valid_splitting = false;
        if(valid_splitting) {
          output.push_back(treeCopy);
        }
      }
    }
  }
  return output;
}

std::vector<Tree<Cluster>> genPermutations(const Tree<Cluster>& tree);

void GenTrees(const Tree<Cluster>& input, int nUnresolved, std::vector<Tree<Cluster>>& output) {
  // Takes an input tree and adds nUresolved particles in possible combinations to the tree. The new trees are stored to the output vector
  if(nUnresolved <= 0) {
    output.push_back(input);
    return;
  }
  int nUnresolved_sofar = 0;
  for(TreeNode<Cluster>* node : input.getNodes()) {
    if((!(node->data.isReference)) and (node != input.getRoot())) nUnresolved_sofar++;
  }
  std::vector<Tree<Cluster>> next_level(GenTrees(input, nUnresolved));
  for(Tree<Cluster> tree : next_level) {
    int nUnresolved_next = 0;
    std::vector<TreeNode<Cluster>*> nodes = tree.getNodes();
    for(TreeNode<Cluster>* node : nodes) {
      if((!(node->data.isReference)) and (node != tree.getRoot())) nUnresolved_next++;
    }
    // Recursively call the function. In this call we generated (nUnresolved_next - nUnresolved_sofar) unresolved nodes,
    // so in the next iteration we only need to generate nUnresolved - (nUnresolved_next - nUnresolved_sofar) unresolved nodes
    GenTrees(tree, nUnresolved - (nUnresolved_next - nUnresolved_sofar), output);
  }
}

void removeRedundentTrees(std::vector<Tree<Cluster>>& trees);

std::vector<Tree<Cluster>> GenTrees(int nUnresolved) {
  // This function generates all possible trees with nUnresolved unresolved nodes.
  Tree<Cluster> baseTree;
  TreeNode<Cluster>* root = new TreeNode<Cluster>(Cluster());
  baseTree.setRoot(root);
  //TreeNode<Cluster>* r1 = new TreeNode<Cluster>(Cluster(true));
  //baseTree.addChild(root, r1);
  std::vector<Tree<Cluster>> output;
  int nReference = 1;
  for(int nReference = 1; nReference <= nUnresolved; nReference++) {
    Tree<Cluster> treeCopy(baseTree);
    for(int rDummy = 1; rDummy <= nReference; rDummy++) {
      TreeNode<Cluster>* r = new TreeNode<Cluster>(Cluster(true));
      treeCopy.addChild(treeCopy.getRoot(), r);
    }
    std::vector<Tree<Cluster>> trees;
    GenTrees(treeCopy, nUnresolved, trees);
    output.insert(output.end(), trees.begin(), trees.end());
  }

  removeRedundentTrees(output);
  // Determine the number of unresolved partons in each cluster
  for(Tree<Cluster>& tree : output) {
    std::vector<TreeNode<Cluster>*> nodes = tree.getNodes();
    for(TreeNode<Cluster>* node : nodes) {
      int nUnresolved = 0;
      for(TreeNode<Cluster>* child : node->children) {
        if(!child->data.isReference) nUnresolved++;
      }
      node->data.unresolved = nUnresolved;
    }
  }
  return output;
}

void compareNodes(TreeNode<Cluster>* node1, TreeNode<Cluster>* node2, bool& status) {
  if(!status) return;
  if(node1 == nullptr and node2 == nullptr) return;
  if(node1->children.size() != node2->children.size()) {
    status = false;
    return;
  }
  else {
    for(int i = 0; i < node1->children.size(); i++) {
      compareNodes(node1->children[i], node2->children[i], status);
    }
  }
}

bool compareTrees(const Tree<Cluster>& tree1, const Tree<Cluster>& tree2) {
  TreeNode<Cluster>* root1 = tree1.getRoot();
  TreeNode<Cluster>* root2 = tree2.getRoot();
  bool output = true;
  if(root1->children.size() != root2->children.size()) {
    return false;
  }
  else {
    for(int i = 0; i < root1->children.size(); i++) {
      if(root1->children[i]->data.isReference != root2->children[i]->data.isReference){
        std::cout << root1->children[i]->data.isReference  << ", " << root2->children[i]->data.isReference << "\t" << root1->children.size() << std::endl;
        return false;
      }
    }
    // If we get here, than the first level must be equal. Now check the remaining levels
    for(int i = 0; i < root1->children.size(); i++) {
      compareNodes(root1->children[i], root2->children[i], output);
    }
  }
  return output;
}

void genPermutations(const Tree<Cluster>& tree, TreeNode<Cluster>* input, std::vector<Tree<Cluster>>& output) {
  // Takes a tree and a node inside that tree, and then generates all possible permutations of the children of the input node.
  // The resulting trees are stored in the output vector
  if(input->children.size() == 0) {
    output.push_back(tree);
    return;
  }
  std::vector<TreeNode<Cluster>*> nodes = tree.getNodes();
  int input_index = -1;
  for(int i = 0; i < nodes.size(); i++) {
    if(nodes[i] == input) input_index = i;
  }
  std::vector<int> range;
  for(int i = 0; i < input->children.size(); i++) range.push_back(i);
  std::vector<std::vector<int>> permutations = getPermutations(range);
  for(std::vector<int>& permutation : permutations) {
    Tree<Cluster> treeCopy(tree);
    TreeNode<Cluster>* inputCopy = treeCopy.getNodes()[input_index];
    std::vector<TreeNode<Cluster>*> nodes_new = inputCopy->children;
    for(int j = 0; j < permutation.size(); j++) {
      inputCopy->children[j] = nodes_new[permutation[j]];
    }
    for(TreeNode<Cluster>* child : inputCopy->children) {
      genPermutations(treeCopy, child, output);
    }
  }
}

std::vector<Tree<Cluster>> genPermutations(const Tree<Cluster>& tree) {
  std::vector<Tree<Cluster>> output;
  genPermutations(tree, tree.getRoot(), output);
  int counter = 1;
  int i = 0;
  while(i < output.size() - 1) {
    int j = i + 1;
    while(j < output.size()) {
      Tree<Cluster> original = output[i];
      Tree<Cluster> permutation = output[j];
      if(compareTrees(original, permutation)) {
        output.erase(output.begin() + j);
      }
      else {
        j++;
      }
    }
    i++;
  }
  return output;
}

void removeRedundentTrees(std::vector<Tree<Cluster>>& trees) {
  int i = 0;
  while(i < trees.size()) {
    Tree<Cluster> tree_i = trees[i];
    std::vector<Tree<Cluster>> permutations = genPermutations(tree_i);
    // Loop through remaining trees and check if one of the trees can be generated through a permutation
    int j = i + 1;
    while(j < trees.size()) {
      Tree<Cluster> tree_j = trees[j];
      // Check if tree_j is in permutations
      bool isPermutation = false;
      for(int k = 1; k < permutations.size(); k++) {
        if(compareTrees(tree_j, permutations[k])) {
          isPermutation = true;
          break;
        }
      }
      if(isPermutation) {
        trees.erase(trees.begin() + j);
      }
      else {
        j++;
      }
    }
    i++;
  }
}

void GenSectors(TreeNode<Cluster>* node, int reference, int& unresolved) {
  // This function labels a single tree according to the indices provided in reference
  if(node == nullptr) return;

  if(node->data.isReference) {
    node->data.reference = reference;
    for(TreeNode<Cluster>* child : node->children) {
      GenSectors(child, reference, unresolved);
    }
  }
  else {
    node->data.reference = unresolved;
    for(TreeNode<Cluster>* child : node->children) {
      unresolved++;
      GenSectors(child, node->data.reference, unresolved);
    }
  }
}

std::vector<Tree<Cluster>> GenSectors(std::vector<int> flavor, const Tree<Cluster>& tree, int nBorn) {
  // This function takes a bare tree and assigns reference indices to the tree, so that it can be used for the generation of phase-space points
  // flavor is a vector of parton flavor: 0 means no QCD interaction, 1 means QCD interaction
  std::vector<Tree<Cluster>> output;
  int nReference = 0;
  for(TreeNode<Cluster>* child : tree.getRoot()->children){
    if(child->data.isReference) nReference++;
  }

  std::vector<int> flavor_indices; // Contains the indices of the QCD interacting particles
  for(int i = 0; i < flavor.size(); i++) {
    if(flavor[i] == 0) {
      continue;
    }
    else {
      flavor_indices.push_back(i);
    }
  }
  Tree<Cluster> treeCopy;
  if(nReference > flavor_indices.size()) return output;
  else if(nReference == flavor_indices.size()) {
    TreeNode<Cluster>* root = new TreeNode<Cluster>(Cluster());
    treeCopy.setRoot(root);
    for(int i = 0; i < nReference; i++) {
      TreeNode<Cluster>* r = new TreeNode<Cluster>(Cluster(true));
      r->data.unresolved = 0;
      treeCopy.addChild(root, r);
    }
    for(int level_counter = 1; level_counter <= nReference; level_counter++) {
      std::vector<TreeNode<Cluster>*> level = treeCopy.getLevel(level_counter);
      TreeNode<Cluster>* u = new TreeNode<Cluster>(Cluster(false));
      u->data.unresolved = 0;
      if(level_counter == 1) {
        level[0]->data.unresolved = 1;
        treeCopy.addChild(level[0], u);
      }
      else {
        level[2]->data.unresolved = 1;
        treeCopy.addChild(level[2], u);
      }
      for(int i = level_counter==1?0:2; i < level.size(); i++) {
        if(level[i]->data.isReference) {
          TreeNode<Cluster>* r = new TreeNode<Cluster>(Cluster(true));
          r->data.unresolved = 0;
          treeCopy.addChild(level[i], r);
        }
      }
    }
  }
  std::vector<std::vector<int>> distributions = generateSubsets(flavor_indices, nReference);
  for(std::vector<int>& distribution : distributions) {
    std::vector<std::vector<int>> permutations = getPermutations(distribution);
    for(std::vector<int> permutation : permutations) {
      Tree<Cluster> treeCopyCopy;
      if(nReference == flavor_indices.size()) {
        treeCopyCopy = Tree<Cluster>(treeCopy);
      }
      else {
        treeCopyCopy = Tree<Cluster>(tree);
      }
      int unresolved_counter = nBorn;
      for(int i = 0; i < treeCopyCopy.getRoot()->children.size(); i++) {
        GenSectors(treeCopyCopy.getRoot()->children[i], permutation[i], unresolved_counter);
      }
      output.push_back(treeCopyCopy);
    }
  }

  return output;
}