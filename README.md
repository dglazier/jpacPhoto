#   jpacPhoto
Framework for amplitude analysis involving single meson production via quasi-elastic scattering of a real photon on a nucleon target. Focus on expandability and easy interfacing with Monte-Carlo tools and event generators.

<p align="center">
  <img width="300" src="./doc/FeynmanDiagram.png">
</p>

Such processes are of interest at many experiments at JLab and the future EIC.

Requires [ROOT](https://root.cern.ch/) (tested with version 6.17) with [*MathMore*](https://root.cern.ch/mathmore-library) libraries installed.

##  INSTALLATION
To install use:
```bash
git clone --recursive https://github.com/dwinney/jpacPhoto.git 
cd jpacPhoto
mkdir build && cd build
cmake ..
cmake --build . --target install
```
This will create a `jpacPhoto/lib/libjpacPhoto.so` and `jpacPhoto/lib/libjpacStyle.so` libraries that can be linked to other code. All executables (e.g. to reproduce plots and results in [[1]](https://arxiv.org/abs/1907.09393) and [[2]](https://arxiv.org/abs/2008.01001)) will also be build in `jpacPhoto/bin`.

Note: cloning with the `--recursive` flag is required to initiate the `jpacStyle` submodule. Compiling can also be sped up with the optional `-jN` flag where N is the number of available cores.

##  AMPLITUDES
The main object of interest is the abstract [`amplitude`](./include/amplitudes/amplitude.hpp) class. This allows you to build [observables](./src/amplitudes/observables.cpp) from helicity amplitudes:

* Probability distribution ( Σ_λ | A |^2 )
* Differential cross section ( dσ / dt )
* Integrated total cross section ( σ )
* Polarization asymmetries ( A_LL and K_LL )
* Spin density matrix elements ( ρ^α_λ,λ' )
* Integrated beam asymmetry ( Σ )
* Parity asymmetry ( P_σ )

Available amplitudes, so far, include:

* [Single baryon resonance](./include/amplitudes/baryon_resonance.hpp) (s-channel)
* [Pomeron exchange](./include/amplitudes/pomeron_exchange.hpp) (t-channel)
* [(fixed-spin and reggeized) Charged pseudo-scalar meson exchange](./include/amplitudes/vector_exchange.hpp) (t-channel)
* [(fixed-spin and reggeized) Vector meson exchange](./include/amplitudes/vector_exchange.hpp) (t-channel)
* [(fixed-spin) Dirac fermion exchange](./include/amplitudes/dirac_exchange.hpp) (u-channel)
* [(fixed-spin) Rarita-Schwinger fermion exchange](./include/amplitudes/rarita_exchange.hpp) (u-channel)

Incoherent (interfering) sums of amplitudes may be constructed through the [`amplitude_sum`](./include/amplitudes/amplitude_sum.hpp) class.

##  EXECUTABLES
The executables folder includes multiple applications of the above amplitudes to different reactions. See documentation in each respective `.cpp` file for usage.

##  PLOTTING
Plots are automatically created using the JPAC collaboration style guidelines. For more information see the [jpacStyle](https://github.com/dwinney/jpacStyle) library.

<p align="center">
  <img width="275" src="./doc/JPAClogo.png">
</p>

##  REFERENCES
+ [1] [Double Polarization Observables in Pentaquark Photoproduction](https://arxiv.org/abs/1907.09393)
+ [2] [XYZ spectroscopy at electron-hadron facilities: Exclusive processes](https://arxiv.org/abs/2008.01001)
+ [3] [JPAC Website](http://cgl.soic.indiana.edu/jpac/index.php)