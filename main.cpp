#include <iostream>
#include "src/PhaseSpace.hpp"
#include "src/Tree.hpp"
#include "/home/tom/Documents/software/software/Cuba-4.2.2/cuba.h"
#include "VEGAS_interface.hpp"


int main() {
  srand(12);
  double COM = 1000.;
  int nBorn = 3;
  int nUnresolved = 0;

  int sample_size = 1000000;
  //cubareal x[3*nBorn - 4];
  int spin = -1;
  int neval;
  int fail;
  cubareal integral[1];
  cubareal error[1];
  cubareal prob[1];
  Vegas(3*nBorn - 4, 1, *integrand_born, &COM, 1, 0.0001, 0.0001, 8, 12, 100, sample_size, 1000, 10000, 1000, 1, "", &spin, &neval, &fail, integral, error, prob);
  std::cout << "Born configuration:" << std::endl;
  std::cout << "integral = " << integral[0] << " +- " << error[0] << "\t" << error[0]/integral[0]*100 << " %\t" << prob[0] << std::endl;
  PhaseSpace control = RAMBO(nBorn, COM);
  std::cout << "integralControl = " << control.weight << "\n\n" << std::endl;

  Cluster cluster1(4, 3);
  Cluster cluster2(4, 2);
  Cluster cluster3(5, 2);
  std::vector<Cluster> cluster = {cluster1};
  for(auto& c : cluster) {
    nUnresolved += c.unresolved;
  }
  TreeNode<Cluster>* rootc = new TreeNode<Cluster>(Cluster(4, 1));
  TreeNode<Cluster>* node1c = new TreeNode<Cluster>(Cluster(4, 1));
  TreeNode<Cluster>* node2c = new TreeNode<Cluster>(Cluster(5, 1));
  Tree<Cluster> clusterTree;
  clusterTree.setRoot(rootc);
  clusterTree.addChild(rootc, node1c);
  clusterTree.addChild(rootc, node2c);
  clusterTree.print();

  UserData data(COM, nBorn, cluster);
  Vegas(3*(nBorn + nUnresolved) - 4, 1, *integrand_full, &data, 1, 0.001, 0.001, 0, 12, 100, sample_size, 1000, 10000, 1000, 2, "", &spin, &neval, &fail, integral, error, prob);
  std::cout << "Full configuration:" << std::endl;
  std::cout << "integral = " << integral[0] << " +- " << error[0] << "\t" << error[0]/integral[0]*100 << " %\t" << prob[0] << std::endl;
  control = RAMBO(nBorn + nUnresolved, COM);
  std::cout << "integralControl = " << control.weight << std::endl;
  std::cout << "ratio = " << integral[0]/control.weight << std::endl;
  std::cout << "deviation = " << (integral[0] - control.weight)/error[0] << std::endl;

  UserData2 data2(COM, nBorn, nUnresolved, &clusterTree);
  Vegas(3*(nBorn + nUnresolved) - 4, 1, *integrand_full2, &data2, 1, 0.001, 0.001, 0, 12, 100, sample_size, 1000, 10000, 1000, 2, "", &spin, &neval, &fail, integral, error, prob);
  std::cout << "\n\nFull configuration2:" << std::endl;
  std::cout << "integral = " << integral[0] << " +- " << error[0] << "\t" << error[0]/integral[0]*100 << " %\t" << prob[0] << std::endl;
  control = RAMBO(nBorn + nUnresolved, COM);
  std::cout << "integralControl = " << control.weight << std::endl;
  std::cout << "ratio = " << integral[0]/control.weight << std::endl;
  std::cout << "deviation = " << (integral[0] - control.weight)/error[0] << std::endl;

  // Test Infrared limits
  PhaseSpace pp = Splitting(nBorn, COM);
  pp.print();
  std::cout << "Check limits ..." << std::endl;
  double scale = 1;
  while (scale > 1.e-8) {
    scale *= 0.1;
    std::vector<std::vector<std::vector<double>>> xParFull;
    for(int c = 0; c < 3; c++) {
      double scale2 = scale;
      if(c != 0) scale2 = scale*scale;
      double eta = scale2;
      double xi = 0.5;
      double phi = rnd(0., 1.);
      std::vector<std::vector<double>> xPar = {{eta, xi, phi}};
      xParFull.push_back(xPar);
    }
    PhaseSpace ppFull = GenMomenta2(pp, clusterTree, xParFull);
    std::cout << scale << ": " << std::endl;
    std::cout << "\t" << ppFull.momenta[4]*ppFull.momenta[5]/std::pow(COM, 2) << "\t" << ppFull.momenta[4]*ppFull.momenta[5]/std::pow(COM*scale, 2) << std::endl;
    std::cout << "\t" << ppFull.momenta[6]*ppFull.momenta[7]/std::pow(COM, 2) << "\t" << ppFull.momenta[6]*ppFull.momenta[7]/std::pow(COM*scale, 2) << std::endl;
    std::cout << "\t" << (ppFull.momenta[4] + ppFull.momenta[5])*(ppFull.momenta[6] + ppFull.momenta[7])/std::pow(COM, 2)
              << "\t" << (ppFull.momenta[4] + ppFull.momenta[5])*(ppFull.momenta[6] + ppFull.momenta[7])/std::pow(COM, 2)/scale << std::endl;

    //ppFull.print();
    std::cout << std::endl;
  }

  // Perform Monte Carlo integration over phase space
  std::cout << "\n################################################################\n" << std::endl;
  int samples = 100000;
  double var = 0., varRAMBO = 0;
  double err = 0., errRAMBO = 0;
  double result = 0., resultRAMBO = 0;
  for(int event_counter = 0; event_counter < samples; event_counter++) {
    //PhaseSpace pp = RAMBO(nBorn, COM);
    PhaseSpace pp = Splitting(nBorn, COM);
    //PhaseSpace event = GenMomenta(pp, cluster);
    PhaseSpace event = GenMomenta2(pp, clusterTree);
    PhaseSpace eventRAMBO = RAMBO(nBorn + nUnresolved, COM);

    /*std::vector<std::vector<double>> xPar;
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
    PhaseSpace eventRAMBO = RAMBO(nBorn, COM);*/
    cubareal f[1];
    cubareal x[3*(nBorn + nUnresolved) - 4];
    for(int i = 0; i < 3*(nBorn + nUnresolved) - 4; i++) {
      x[i] = rnd(0, 1);
    }
    int ncomp = 1;
    const int ndim = 3*(nBorn + nUnresolved) - 4;
    //int dummy = integrand_full(&ndim, x, &ncomp, f, &data);
    //event.weight = f[0];


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