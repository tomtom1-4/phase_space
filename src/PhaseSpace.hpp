#ifndef PHASESPACE_HPP
#define PHASESPACE_HPP

#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include "Exception.hpp"
#include "Utilities.hpp"

class LorentzMatrix {
  public:
    std::vector<std::vector<double>> components = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    LorentzMatrix(std::vector<std::vector<double>> components) : components(components){}
    LorentzMatrix(){};
    void print();
    void print_list();
};

class Momentum {
  public:
    std::vector<double> components = {0,0,0,0};

    void print() const {
      std::cout << std::setprecision(15)  << "{" << components[0] << ", " << components[1] << ", " << components[2] << ", " << components[3] << "}" << std::endl;
    }

    double inv_mass() const {
      try {
        return std::sqrt(components[0]*components[0] - components[1]*components[1] - components[2]*components[2] - components[3]*components[3]);
      }
      catch(const Exception& e) {
        std::cerr << e.what() << " Could not evaluate invariant mass of momentum." << std::endl;
      }
    }

    double inv_mass2() const {
      try {
        return components[0]*components[0] - components[1]*components[1] - components[2]*components[2] - components[3]*components[3];
      }
      catch(const Exception& e) {
        std::cerr << e.what() << " Could not evaluate invariant mass of momentum." << std::endl;
      }
    }

    Momentum(std::vector<double> components) : components(components){}

    Momentum(){};
};

Momentum operator+(Momentum p1, Momentum p2);

Momentum operator-(Momentum p1, Momentum p2);

double operator*(Momentum p1, Momentum p2);

Momentum operator*(double a, Momentum p);

Momentum operator*(Momentum p, double a);

Momentum operator/(Momentum p, double a);

Momentum operator*(LorentzMatrix lam, Momentum p);

Momentum operator-(Momentum p);

LorentzMatrix operator*(LorentzMatrix lam1, LorentzMatrix lam2);

class PhaseSpace {
  public:
    std::vector<Momentum> momenta;

    PhaseSpace(std::vector<Momentum> momenta) : momenta(momenta){}

    PhaseSpace(){}

    int size() const {
      return momenta.size();
    }

    void check_momentum_conservation(double acc = 1.e-8);

    void check_onshellness(double acc);

    void print();

    double weight;

};

class Cluster {
  public:
    int reference; // index of the reference momentum
    int unresolved; // number of unresolved partons
    Cluster(int reference, int unresolved) : reference(reference), unresolved(unresolved) {};
    std::vector<Momentum> unresolved_momenta;
    Momentum reference_momentum;
    double weight;
};

class PESCPhaseSpace : public PhaseSpace  {
  public:
    std::vector<Cluster> cluster;
};

double rnd(double lower, double upper);

double RAMBO_measure(int nMomenta, double COM);

PhaseSpace RAMBO(const int nMomenta, const double COM);

PhaseSpace Splitting(int nMomenta, double COM, std::vector<std::vector<double>> x);

PhaseSpace Splitting(int nMomenta, double COM);

PESCPhaseSpace GenMomenta(const PhaseSpace pp, const std::vector<Cluster>& cluster, std::vector<std::vector<std::vector<double>>> x);

PESCPhaseSpace GenMomenta(const PhaseSpace pp, const std::vector<Cluster>& cluster);

#endif