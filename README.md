README 

This package provides the software described in  [Aussie Greene and Larry L. Schumaker, 
Solving PDEs on Curved Domains in 3D using the Immersed Penalized Boundary Method 
with Trivariate Spline Spaces] for  numerically solving second order elliptic
boundary-value problems defined on 3D domains with curved boundaries. The package
contains the Matlab implementations of the functions and scripts described
in the paper, along with certain supporting functions for computing with splines
that are described in the book [Larry L. Schumaker, Spline Functions: Computational
Methods, SIAM, 2015].

Each example in the paper can be reproduced by running the designated script
given in the example. At run time the user will be prompted for various inputs which 
are also listed in the example. There are 17 scripts whose names begin with T3.
They work with trivariate polynomial splines on tetrahedral partitions. There are
an additional 13 scripts whose names begin with TP3 -- they work with 
trivariate tensor-product splines. All of the scripts produce one or more figures. 

The package also contains several files with the extensions .dat, .stl, and
.ply.  These are used to input the domains discussed in the paper. To read in
one of this files the user simply responds to the first prompt in each script with 
the name of file to be read. The functions setpde3.m and  setpde3L.m are 
used to set up the boundary-value problem. They are called with a parameter pde 
which selects the desired problem. The second prompt in each script asks the user
for this parameter.
