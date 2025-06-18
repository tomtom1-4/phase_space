#ifndef PHASESPACE2_HPP
#define PHASESPACE2_HPP

#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <variant>
#include "Exception.hpp"
#include "Utilities.hpp"
#include "Tree.hpp"

namespace PSF
{
/**
 * @class LorentzMatrix
 * @brief Represents a 4×4 matrix commonly used in Lorentz transformations.
 *
 * This class provides a simple 4×4 matrix implementation with double-precision
 * floating point entries, intended for use in special relativity calculations.
 * By default, it is initialized to a zero matrix.
 *
 * Example usage:
 * @code{.cpp}
 * LorentzMatrix mat({
 *     {1.0, 0.0, 0.0, 0.0},
 *     {0.0, 1.0, 0.0, 0.0},
 *     {0.0, 0.0, 1.0, 0.0},
 *     {0.0, 0.0, 0.0, 1.0}
 * });
 * mat.print();
 * @endcode
 */
class LorentzMatrix {
public:
  /**
   * @brief Stores the 4×4 components of this matrix.
   *
   * Each row is a std::vector of 4 doubles. By default, the matrix is
   * initialized to all zeros.
   */
  std::vector<std::vector<double>> components = {
      {0, 0, 0, 0},
      {0, 0, 0, 0},
      {0, 0, 0, 0},
      {0, 0, 0, 0}
  };

  /**
   * @brief Constructs a LorentzMatrix from a 4×4 arrangement of doubles.
   * @param components A 4×4 vector of doubles to initialize the matrix with.
   */
  LorentzMatrix(std::vector<std::vector<double>> components)
      : components(components) {}

  /**
   * @brief Default constructor initializing the matrix to zero.
   */
  LorentzMatrix() = default;

  /**
   * @brief Prints the matrix in a neatly formatted 4×4 grid.
   */
  void print() const;

  /**
   * @brief Prints the matrix components in a single list format.
   */
  void print_list() const;
};

/**
 * @class Momentum
 * @brief Represents a 4-momentum in special relativity.
 *
 * The four components are stored in the order {E, px, py, pz}.
 * By default, they are initialized to zero.
 * This class provides methods for printing the current momentum
 * and for calculating its invariant mass and invariant mass squared.
 */
class Momentum {
public:
  /**
   * @brief The four components of the 4-momentum.
   *
   * Stored as {E, px, py, pz}.
   * Defaults to {0, 0, 0, 0}.
   */
  std::vector<double> components = {0, 0, 0, 0};

  /**
   * @brief Prints the 4-momentum components to standard output.
   *
   * Displays high precision using std::setprecision(15).
   * The format is {E, px, py, pz}.
   */
  void print() const {
      std::cout << std::setprecision(15)
                << "{" << components[0] << ", "
                << components[1] << ", "
                << components[2] << ", "
                << components[3] << "}" << std::endl;
  }

  /**
   * @brief Calculates and returns the invariant mass.
   *
   * Uses the metric signature (+, -, -, -): \f$m = \sqrt{E^2 - p_x^2 - p_y^2 - p_z^2}\f$.
   *
   * @return Invariant mass (double). If the quantity inside the square root is negative,
   *         the behavior depends on the standard library's sqrt function (NaN for negative input).
   * @note This function throws any exceptions that may arise from the std::sqrt call
   *       (e.g., domain error) if the negative argument is encountered;
   *       the catch block logs the error message.
   */
  double inv_mass() const {
    try {
      return std::sqrt(components[0]*components[0]
                      - components[1]*components[1]
                      - components[2]*components[2]
                      - components[3]*components[3]);
    }
    catch(const Exception& e) {
      std::cerr << e.what() << " Could not evaluate invariant mass of momentum." << std::endl;
    }
  }

  /**
   * @brief Calculates and returns the invariant mass squared.
   *
   * Uses the metric signature (+, -, -, -): \f$m^2 = E^2 - p_x^2 - p_y^2 - p_z^2\f$.
   *
   * @return Invariant mass squared (double). It can be negative if the momentum is spacelike.
   * @note This function also uses a try-catch block. If an exception is thrown,
   *       it will be logged.
   */
  double inv_mass2() const {
    try {
      return components[0]*components[0]
            - components[1]*components[1]
            - components[2]*components[2]
            - components[3]*components[3];
    }
    catch(const Exception& e) {
      std::cerr << e.what() << " Could not evaluate invariant mass of momentum." << std::endl;
    }
  }

  /**
   * @brief Constructs a 4-momentum with specified components.
   * @param components A vector of 4 doubles in the order {E, px, py, pz}.
   */
  Momentum(std::vector<double> components) : components(components) {}

  /**
   * @brief Default constructor initializing the 4 components to zero.
   */
  Momentum() = default;
};

/**
  * @brief Adds two Momentum objects component-wise.
  *
  * @param p1 The first Momentum object.
  * @param p2 The second Momentum object.
  * @return A new Momentum object containing the sum of p1 and p2.
  */
Momentum operator+(Momentum p1, Momentum p2);

/**
  * @brief Subtracts one Momentum from another component-wise.
  *
  * @param p1 The first Momentum object.
  * @param p2 The second Momentum object.
  * @return A new Momentum object containing the difference p1 - p2.
  */
Momentum operator-(Momentum p1, Momentum p2);

/**
  * @brief Computes the Minkowski (dot) product of two 4-momenta.
  *
  * Uses the metric signature (+, - , -, -).
  * @param p1 The first Momentum object.
  * @param p2 The second Momentum object.
  * @return The scalar product (double).
  */
double operator*(Momentum p1, Momentum p2);

/**
  * @brief Multiplies a Momentum by a scalar from the left side (a * p).
  *
  * @param a The scalar factor.
  * @param p The Momentum object.
  * @return A new Momentum object with components scaled by a.
  */
Momentum operator*(double a, Momentum p);

/**
  * @brief Multiplies a Momentum by a scalar from the right side (p * a).
  *
  * @param p The Momentum object.
  * @param a The scalar factor.
  * @return A new Momentum object with components scaled by a.
  */
Momentum operator*(Momentum p, double a);

/**
  * @brief Divides a Momentum by a scalar (p / a).
  *
  * @param p The Momentum object.
  * @param a The scalar divisor.
  * @return A new Momentum object with each component divided by a.
  */
Momentum operator/(Momentum p, double a);

/**
  * @brief Transforms a Momentum by a LorentzMatrix (lam * p).
  *
  * @param lam The LorentzMatrix representing the transformation.
  * @param p The Momentum object.
  * @return A new Momentum object that is the result of applying lam to p.
  */
Momentum operator*(LorentzMatrix lam, Momentum p);

/**
  * @brief Returns the unary negation of a Momentum (i.e., flips the sign of each component).
  *
  * @param p The Momentum object to be negated.
  * @return A new Momentum object with negated components.
  */
Momentum operator-(Momentum p);

/**
  * @brief Multiplies two LorentzMatrix objects (lam1 * lam2).
  *
  * Represents matrix multiplication of 4×4 matrices.
  * @param lam1 The left-hand LorentzMatrix.
  * @param lam2 The right-hand LorentzMatrix.
  * @return A new LorentzMatrix that is the result of lam1 × lam2.
  */
LorentzMatrix operator*(LorentzMatrix lam1, LorentzMatrix lam2);

/**
  * @class PhaseSpace
  * @brief Represents a collection of momenta in phase space, along with utilities for consistency checks.
  *
  * This class stores a vector of Momentum objects and provides functions
  * to check momentum conservation, on-shell conditions, and to print state information.
  */
class PhaseSpace {
public:
  /**
    * @brief Holds a list of Momentum objects.
    *
    * Represents particles or other physical entities in the phase space.
    */
  std::vector<Momentum> momenta;

  /**
    * @brief Constructs a PhaseSpace with the specified list of momenta.
    *
    * @param momenta A vector of Momentum objects to store.
    */
  PhaseSpace(std::vector<Momentum> momenta) : momenta(momenta){}

  /**
    * @brief Default constructor initializing an empty phase space.
    */
  PhaseSpace(){}

  /**
    * @brief Returns the number of Momentum objects in this phase space.
    *
    * @return The size of the momenta vector.
    */
  int size() const {
    return momenta.size();
  }

  /**
    * @brief Checks total momentum conservation within a specified tolerance.
    *
    * Ensures that the sum of all 4-momenta is zero (or close to zero) within the given threshold.
    * @param acc The allowed numerical tolerance (default is 1.e-8).
    * @throws Any exception or warning mechanism if momentum is not conserved within acc.
    */
  void check_momentum_conservation(double acc = 1.e-8);

  /**
    * @brief Checks if each momentum is on-shell within a specified tolerance.
    *
    * Uses the condition \f$ E^2 - p^2 = m^2 \f$ or a similar approach
    * depending on your definition of on-shell.
    * @param acc The allowed numerical tolerance (default is 1.e-5).
    * @throws Any exception or warning mechanism if one or more momenta are off-shell beyond acc.
    */
  void check_onshellness(double acc = 1.e-5);

  /**
    * @brief Prints the phase space information.
    *
    * Includes printing each Momentum in momenta as well as other details such as weight if desired.
    */
  void print() const;

  /**
    * @brief Holds a weight factor associated with this phase space configuration.
    *
    * This can be used in Monte Carlo integrations or event generation.
    */
  double weight;
};

class Cluster {
  public:
    int reference; // index of the reference momentum
    int unresolved; // number of unresolved partons
    bool isReference;
    Cluster(int reference, int unresolved) : reference(reference), unresolved(unresolved) {};
    Cluster(bool isReference) : isReference(isReference) {};
    Cluster() {};
    std::vector<Momentum> unresolved_momenta;
    Momentum reference_momentum;

    bool operator < (const Cluster& cluster2) const {
      return (this->unresolved < cluster2.unresolved);
    }
};

class PESCPhaseSpace : public PhaseSpace  {
  public:
    std::vector<Cluster> cluster;
};

/**
  * @brief Returns random double in [lower, upper] uniquely distributed
  *
  * @param lower lower bound
  * @param upper upper bound
  */
double rnd(double lower, double upper);

/**
  * @brief Returns phase-space point using the RAMBO algorithm
  *
  * The phase space is generated with perfect weight, i.e. every phase-space point has the same weight
  * Kleiss, R., Stirling, W., & Ellis, S. (1986). A new Monte Carlo treatment of multiparticle phase space at high energies. Computer Physics Communications, 40(2–3), 359–373.
  * @param nMomenta Number of final-state momenta to generate.
  * @param COM Center of mass energy.
  */
PhaseSpace RAMBO(const int nMomenta, const double COM);

/**
 * @brief Returns phase-space point using a recursive Splitting algorithm.
 *
 * In every recursion step, the phase space is split into a massless external momentum
 * and an off-shell momentum of mass x[i]*M[i-1], where M[i-1] is the invariant mass of
 * momentum configuration at the previous step.
 *
 * \verbatim
 *       /      /      /      /
 *  COM /  M1  /  M2  /  M3  /
 * ==========================-----
 * \endverbatim
 *
 * @param nMomenta Number of final-state momenta to generate.
 * @param COM Center of mass energy.
 * @param x Random variables used to generate the momenta.
 *          Structure is {{xM_1, xM_2,...,xM_{nMomenta - 2}},
 *                        {xCos_1, xCos_2,...,xCos_{nMomenta - 1}},
 *                        {xPhi_1, xPhi_2,...,xPhi_{nMomenta - 2}}}
 */
PhaseSpace Splitting(int nMomenta, double COM, std::vector<std::vector<double>> x);

PhaseSpace Splitting(int nMomenta, double COM);

PESCPhaseSpace GenMomenta(const PhaseSpace pp, const std::vector<Cluster>& cluster, std::vector<std::vector<std::vector<double>>> x);

PESCPhaseSpace GenMomenta(const PhaseSpace pp, const std::vector<Cluster>& cluster);

PhaseSpace GenMomenta2(const PhaseSpace pp, const Tree<Cluster>& clusterTree, std::vector<std::vector<std::vector<double>>> x);

PhaseSpace GenMomenta2(const PhaseSpace pp, const Tree<Cluster>& clusterTree);

std::vector<Tree<Cluster>> GenTrees(int nUnresolved);

std::vector<Tree<Cluster>> GenSectors(std::vector<int> flavor, const Tree<Cluster>& tree, int nBorn);

bool compareTrees(const Tree<Cluster>& tree1, const Tree<Cluster>& tree2);

}
#endif
