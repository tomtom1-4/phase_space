# Phase Space

Phase Space uses a generalization of the algorithm presented in `[1]` to generate phase-space points. It can generate phase-space points with traditional algorithms like RAMBO, and then enables the user to generate additional unresolved momenta on top of the existing ones parametrized by infrared sensitive variables. It is the first step towards an infrared subtraction scheme at all orders.

## Installation

1. Clone the repository:
```bash
git clone git@github.com:tomtom1-4/phase_space.git
cd phase_space
```

2. Configure: If you want to use phase space integrations with the VEGAS algorithm you can link the project to [CUBA](https://feynarts.de/cuba/). In that case, you have to add the environment variable
```bash
export CUBA_PATH=/path/to/your/CUBA/installation/
```
If you do not want to include CUBA, change ```USE_CUBA``` to false in the ```Makefile```.

3. Build the library:
```bash
make
```
This will build a shared object file as well as the examples in the ```example``` directory.


## Usage

To include the library in any project, ensure that the runtime linker can find the library, e.g. by setting the ```LD_LIBRARY_PATH```
```bash
export LD_LIBRARY_PATH=/path/to/phase_space:$LD_LIBRARY_PATH
```
Include the header files you need in your project
```c++
#include "/path/to/phase_space/src/Utilities.hpp"
#include "/path/to/phase_space/src/Tree.hpp"
#include "/path/to/phase_space/src/PhaseSpace.hpp"
```
Then during compilation, simply include the path to the library
```bash
g++ -o your_project your_project.cpp -L/path/to/phase_space -lphase_space
```

## Getting started
There is a short manual describing the basic algorithm for the generation of phase-space points in the ```docs``` directory.

In the ```examples``` subdirectory you find two example source files. To compile them, first source the environment
```bash
source env.sh
```
then build the target files either individually
```bash
make pp_integration.exe
make gen_pp.exe
```
or all at once
```bash
make
```
If you don't want to include [CUBA](https://feynarts.de/cuba/), then change ```USE_CUBA``` to false in the Makefile of the examples directory.
You may now execute the two example files.

### gen_pp

```gen_pp``` demonstrates how to generate phase-space points. ```phase_space``` has implemented multiple options for the generation of phase-space points. Phase-space points with resolved momenta (Born phase-space point) can be generated with:
1. the RAMBO algorithm
```c++
PhaseSpace RAMBO(int nMomenta, double COM);
```
2. a splitting algorithm
```c++
PhaseSpace Splitting(int nMomenta, double COM);
PhaseSpace Splitting(int nMomenta, double COM, std::vector<std::vector<double>> x);
```

Once we generated a Born phase-space point, we can add an arbitrary number of unresolved momenta parametrized by infrared variables ($\eta, \xi, \phi$). These variables are defined with respect to reference momenta. To specify against which momenta they are defined against, we use a tree like structure of clusters. These trees can be generated manually and automatically. With
```c++
   std::vector<PSF::tree<PSF::Cluster>>> PSF::GenTrees(int nUnresolved);
```
you can generate all possible trees for a given number of unresolved momenta. ```phase_space``` also allows you to interpret these trees graphically through
```c++
void PSF::tree<PSF::Cluster>>::print();
```
At this point the trees are ``unpopulated``, i.e. they have no knowledge of the underlying Born phase-space point. To make the connection, you can assign momentum indices via
```c++
std::vector<PSF::Tree<PSF::Cluster>> PSF::GenSectors(std::vector<bool> flavor, PSF::Tree<PSF::Cluster> tree, int nBorn);
```
Finally, once a sector has been selected, i.e. a tree has been populated, we can generate a phase-space point with
```c++
PSF::PhaseSpace PSF::GenMomenta(PhaseSpace pp, PSF::Tree<PSF::Cluster> tree, std::vector<std::vector<std::vector<double>>> xPar); // Specific parametrization
PSF::PhaseSpace PSF::GenMomenta(PhaseSpace pp, PSF::Tree<PSF::Cluster> tree); // Random parametrization
```

You can take a look at the resulting phase-space points with
```c++
void PSF::PhaseSpace::print() const;
```

### pp_integration

Once you generated a phase-space point, ```phase_space``` will have already computed the corresponding weight in the background. ```pp_integration``` demonstrates how one may use these weights to perform Monte-Carlo integrations. In this example, we simply integrate over the identity, i.e. we compute the phase-space volume, as it is easily calculable analytically and thus allows us to verify our integration.

```phase_space``` also allows performing the Monte-Carlo integration with important sampling using the VEGAS algorithm through an interface to [CUBA](https://feynarts.de/cuba/). Examples for such an interface can be found in ```VEGAS_interface```.

### phase_space

In the ```tests``` subdirectory, you can also find the source file ```phase_space.cpp```. You can compile it in the same manner as the example files. The program tests if the phase-space point is generated successfully, it verifies if the weights are computed correctly, by computing the phase-space volume with a Monte-Carlo integration, and it also tries to parametrize certain infrared limits.

## Doxygen

```phase_space``` uses [Doxygen](https://www.doxygen.nl/manual/install.html) documentation. To get the documentation, install [Doxygen](https://www.doxygen.nl/manual/install.html), then go to the ```docs``` subdirectory and run
```bash
doxygen Doxyfile
```
This will generate documentation for the library in a PDF (```latex/refman.pdf```) and an HTML format (```index.html```).

## References

- `[1]` Czakon, M., Van Hameren, A., Mitov, A., & Poncelet, R. (2019). Single-jet inclusive rates with exact color at $\mathcal{O} ( {\alpha}_s^4 )$. Journal of High Energy Physics, 2019(10). https://doi.org/10.1007/jhep10(2019)262