Parameters for the parametrized EOS are given in the file 
`peos_parameters.dat'  

---
            3  8.47187d+14     : nphase, rhoini_cgs
  5.01187d+14  2.42103d+34     : rho_0, pre_0
  1.00000d+18  2.85100d+00     : rho, gamma
  1.00000d+15  2.98800d+00     : rho, gamma
  5.01187d+14  3.00500d+00     : rho, gamma
  1.00000d+00  3.00500d+00     : rho, gamma
---

nphase : the number of polytropic segments.  
rhoini_cgs   : the rest mass density at the center of a NS, which used 
               only in the TOV solver for a initial guess.  
rho_0, pre_0 : values of the rest mass density and corresponding pressure 
               used to specify the scale of EOS (an adiabatic constant).  
rho    : dividing densities between polytropic segments.
gamma  : value of adiabadic index in the segment, whose 
         upper bound is the corresponding rho.  
