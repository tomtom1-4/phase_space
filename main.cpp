#include <iostream>
#include "src/PhaseSpace.hpp"
#include "/home/tom/Documents/software/software/Cuba-4.2.2/cuba.h"

/*int cuba_example(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata) {
  f[0] = std::exp(-x[0]*x[0])/std::sqrt(M_PI);
  //std::cout << x[0] << " => " << f[0] << std::endl;
  return 0;
}*/

int integrand(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata) {
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


int main() {
  srand(12);
  double COM = 1000.;
  int nBorn = 5;
  int nUnresolved = 0;
  std::vector<std::vector<double>> xPar;
  {
  std::vector<double> M, cos, phi;
  for(int i = 0; i < nBorn - 2; i++) M.push_back(rnd(0., 1.));
  for(int i = 0; i < nBorn - 1; i++) cos.push_back(rnd(0., 1.));
  for(int i = 0; i < nBorn - 1; i++) phi.push_back(rnd(0., 1.));
  xPar.push_back(M);
  xPar.push_back(cos);
  xPar.push_back(phi);
  }
  PhaseSpace pp = Splitting(nBorn, COM, xPar);

  //PhaseSpace pp = RAMBO(nBorn, COM);
  std::cout << "Born Phase-Space Point:" << std::endl;
  pp.print();
  std::cout << "\nCheck Momentum Conservation:" << std::endl;
  pp.check_momentum_conservation(1.e-8);
  std::cout << "\nCheck On-Shellness:" << std::endl;
  pp.check_onshellness(1.e-8);

  int sample_size = 100000;
  cubareal x[3*nBorn - 4];
  int spin = -1;
  int neval;
  int fail;
  cubareal integral[1];
  cubareal error[1];
  cubareal prob[1];
  Vegas(3*nBorn - 4, 1, *integrand, &COM, 1, 0.00001, 0.00001, 8, 12, 100, sample_size, 1000, 1000, 1000, 1, "", &spin, &neval, &fail, integral, error, prob);
  std::cout << "integral = " << integral[0] << " +- " << error[0] << std::endl;
  PhaseSpace control = RAMBO(nBorn, COM);
  std::cout << "integralControl = " << control.weight << std::endl;

  Cluster cluster1(2, 1);
  Cluster cluster2(3, 2);
  Cluster cluster3(4, 1);
  Cluster cluster4(5, 3);
  std::vector<Cluster> cluster = {cluster1};
  for(auto& c : cluster) {
    nUnresolved += c.unresolved;
  }

  PESCPhaseSpace pp_new = GenMomenta(pp, cluster);
  std::cout << "\nFull Phase-Space Point" << std::endl;
  pp_new.print();
  std::cout << "\nCheck Momentum Conservation:" << std::endl;
  pp_new.check_momentum_conservation(1.e-8);
  std::cout << "\nCheck On-Shellness:" << std::endl;
  pp_new.check_onshellness(1.e-8);

  srand(10);
  pp_new = GenMomenta(pp, cluster);
  std::cout << "\nFull Phase-Space Point" << std::endl;
  pp_new.print();
  std::cout << "\nCheck Momentum Conservation:" << std::endl;
  pp_new.check_momentum_conservation(1.e-8);
  std::cout << "\nCheck On-Shellness:" << std::endl;
  pp_new.check_onshellness(1.e-8);

  // Perform Monte Carlo integration over phase space
  std::cout << "\n################################################################\n" << std::endl;
  //double pp_Volume = RAMBO_measure(nBorn + nUnresolved, COM);
  //pp_Volume = 1./4./std::pow(2.*M_PI, 1) * std::pow(COM/2., 2);




  int samples = 100000;
  double var = 0., varRAMBO = 0;
  double err = 0., errRAMBO = 0;
  double result = 0., resultRAMBO = 0;
  for(int event_counter = 0; event_counter < samples; event_counter++) {
    //PhaseSpace pp = RAMBO(nBorn, COM);
    //PhaseSpace event = GenMomenta(pp, cluster);
    //PhaseSpace eventRAMBO = RAMBO(nBorn + nUnresolved, COM);
    std::vector<std::vector<double>> xPar;
    {
    std::vector<double> M, cos, phi;
    for(int i = 0; i < nBorn - 2; i++) M.push_back(rnd(0., 1.));
    for(int i = 0; i < nBorn - 1; i++) cos.push_back(rnd(0., 1.));
    for(int i = 0; i < nBorn - 1; i++) phi.push_back(rnd(0., 1.));
    xPar.push_back(M);
    xPar.push_back(cos);
    xPar.push_back(phi);
    }

    PhaseSpace event = Splitting(nBorn, COM, xPar);
    PhaseSpace eventRAMBO = RAMBO(nBorn, COM);

    double integrand = 1;//std::exp(-std::pow(event.momenta[3]*event.momenta[5], 2));
    double integrandRAMBO = 1;//std::exp(-std::pow(eventRAMBO.momenta[3]*eventRAMBO.momenta[5], 2));


    result = result*event_counter/(event_counter + 1.) + event.weight*integrand/(event_counter + 1.);
    var = var*event_counter/(event_counter + 1.) + std::pow(event.weight*integrand - result, 2)/(event_counter + 1.);

    err = std::sqrt(var)/std::sqrt(event_counter);

    resultRAMBO = resultRAMBO*event_counter/(event_counter + 1.) + eventRAMBO.weight*integrandRAMBO/(event_counter + 1.);
    varRAMBO = varRAMBO*event_counter/(event_counter + 1.) + std::pow(eventRAMBO.weight*integrandRAMBO - resultRAMBO, 2)/(event_counter + 1.);

    errRAMBO = std::sqrt(varRAMBO)/std::sqrt(event_counter);


    if((event_counter + 1) % (samples/10) == 0) {
      std::cout << std::setprecision(5) << event_counter + 1 << "/" << samples << ": " << result << " +- " << err << "\t" << err/result*100
                                                << "%\tRAMBO: " << resultRAMBO << " +- " << errRAMBO << "\t" << errRAMBO/resultRAMBO*100
      << "%\t deviation: " << (result - resultRAMBO)/sqrt(std::pow(err, 2) + std::pow(errRAMBO, 2)) << "\tratio: " << result/resultRAMBO << std::endl;
    }
  }
  std::cout << "\n\ncontrol:" << std::endl;
  std::cout << resultRAMBO << std::endl;
  std::cout << "deviation = " << (result - resultRAMBO)/sqrt(std::pow(err, 2) + std::pow(errRAMBO, 2)) << std::endl;
  std::cout << "ratio = " << result/resultRAMBO << std::endl;
}