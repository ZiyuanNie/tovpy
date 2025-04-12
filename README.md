# Tolman-Oppenheimer-Volkoff solver in Python

This repository provides a package that solves the
Tolman–Oppenheimer–Volkoff (TOV) equations for general-relativistic
nonrotating equilibrium star configurations together with even/odd
stationary perturbations of arbitrary multipolar index. 

Note:
- TOV equations are formulated following [Lindblom
(1992)](https://articles.adsabs.harvard.edu/pdf/1992ApJ...398..569L) 
- Perturbations equations and Love number computation follow [Damour
and Nagar (2009)](https://doi.org/10.1103/PhysRevD.80.084035) 
- The code supports EOS in piecewise polytropic and tabular form

See [this file](AUTHORS.md) for authors and contributors.

The code is distributed under the [GNU GPLv3
licence](https://www.gnu.org/licenses/gpl-3.0.en.html). 

TOVpy has benefited from open source code and libraries, in particular:

 * [MATLAB TOVL code](https://bitbucket.org/bernuzzi/tov) 
 * [LIGO-LAL](https://lscsoft.docs.ligo.org/lalsuite/lalsimulation)
 * [bilby](https://lscsoft.docs.ligo.org/bilby)

