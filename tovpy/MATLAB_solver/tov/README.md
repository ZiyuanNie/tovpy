# Matlab code for general relativistic spherical star solutions

This code solves TOV equations in general relativity together with
even/odd stationary perturbations for barotropic fluids and for
aribitrary multipolar index ell.

The code has been developed from 2007 based on a early version of
Alessandro Nagar. The current version is essentially the 2009
implementation for computing the Love numbers.

Copyright (c) 2009 S.Bernuzzi, A.Nagar

The code is provided "as is", without warranty of any kind.
We ask you to cite the following two papers if you use this code or
any part of it.

```
@article{Damour:2009vw,
      author         = "Damour, Thibault and Nagar, Alessandro",
      title          = "{Relativistic tidal properties of neutron stars}",
      journal        = "Phys. Rev.",
      volume         = "D80",
      year           = "2009",
      pages          = "084035",
      doi            = "10.1103/PhysRevD.80.084035",
      eprint         = "0906.0096",
      archivePrefix  = "arXiv",
      primaryClass   = "gr-qc",
      SLACcitation   = "%%CITATION = ARXIV:0906.0096;%%"
}
```

```
@article{Bernuzzi:2008fu,
      author         = "Bernuzzi, Sebastiano and Nagar, Alessandro",
      title          = "{Gravitational waves from pulsations of neutron stars
                        described by realistic Equations of State}",
      journal        = "Phys. Rev.",
      volume         = "D78",
      year           = "2008",
      pages          = "024024",
      doi            = "10.1103/PhysRevD.78.024024",
      eprint         = "0803.3804",
      archivePrefix  = "arXiv",
      primaryClass   = "gr-qc",
      SLACcitation   = "%%CITATION = ARXIV:0803.3804;%%"
}
```

Note the convention: matlab functions filenames start with Uppercase,
lowercase is used for scripts. 