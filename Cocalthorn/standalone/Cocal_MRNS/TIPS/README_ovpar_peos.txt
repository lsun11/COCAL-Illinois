`ovpar_peos.dat' is used for a computation of TOV solution.
TOV solver is in the `Cocal/spherical_data/TOV_schwartzshild_coord' 
directory, which solves TOV equations using 4th order Runge-Kutta method.
The parameter file `ovpar_peos.dat' is as follows:
---
  400   16  8.3102833E+00        : nstep, ndiv, radiini
    0    1  0.1700000d-00        : itype, iter_mode, chope

itype = 0  iteration for compactness, or mass
itype = 1  single structure
iter_mode = 0  iteration for compactness
iter_mode = 1  iteration for rest mass
iter_mode = 2  iteration for gravitational mass
---

nstep : number of grid points from the center to the surface of NS, where 
        the data of density etc. is returned from the Runge-Kutta routine.
ndiv  : Runge-Kutta routine further divide ndiv points between the grid. 
radini: Initial guess for the radius of NS in [km].
chope is not used.
