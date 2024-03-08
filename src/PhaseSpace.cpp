#include "PhaseSpace.hpp"

void LorentzMatrix::print() {
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

void LorentzMatrix::print_list() {
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

void PhaseSpace::print() {
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
  if ( std::abs(m1 - m2) > 1.e-5) {
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


PESCPhaseSpace GenMomenta(const PhaseSpace pp, const std::vector<Cluster>& cluster) {
  PESCPhaseSpace pp_new;
  pp_new.cluster = cluster;
  std::vector<double> x(pp_new.cluster.size(), 0.);
  int dim = 4;
  double jacobian = 1.;
  Momentum P = (-1.)*(pp.momenta[0] + pp.momenta[1]);
  Momentum rTot, rWeighted, uTot;
  for(int j = 0; j < pp_new.cluster.size(); j++) {
    Momentum r = pp.momenta[pp_new.cluster[j].reference];
    // Find rotation that yields r from (1, 0, 0, 1) and apply to uhat
    LorentzMatrix rotation = find_LT(r/r.components[0], Momentum(std::vector<double>({1., 0., 0., 1.})));
    rTot = rTot + r;
    // Generate the unresolved momenta in the cluster
    for(int a = 0; a < pp_new.cluster[j].unresolved; a++) {
      double eta = rnd(0., 1.);
      double xi = rnd(0., 1.);
      double phi = rnd(0., 1.)*2.*M_PI;
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
      uTot = uTot + u;
    }
    double xj = (-2.*P*rTot + rTot*rTot - uTot*uTot + 2.*P*uTot + 2.*P*rWeighted - rWeighted*rWeighted - 2.*rWeighted*uTot)/
        (-2.*P*r + 2.*r*rWeighted + 2.*r*uTot);
    x[j] = xj;

    // determine weight factor
    double c1 = -2.*P*r + 2.*r*rTot;
    double c2 = -2.*P*r + 2.*r*rWeighted + 2.*r*uTot;
    double c3 = -2.*P*rWeighted + rWeighted*rWeighted - 2.*P*uTot + uTot*uTot + 2.*rWeighted*uTot + 2.*P*(rTot - r) + (rTot - r)*(rTot - r);

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

  pp_new.momenta.push_back(pp.momenta[0]);
  pp_new.momenta.push_back(pp.momenta[1]);

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
      pp_new.momenta.push_back(pi_new);
    }
    else {
      pp_new.momenta.push_back(pp_new.cluster[cluster_index].reference_momentum);
      for(auto& u:pp_new.cluster[cluster_index].unresolved_momenta) {
        pp_new.momenta.push_back(u);
      }
    }
  }
  pp_new.weight = pp.weight*jacobian;
  return pp_new;
}