In the following, contents of the parameter file 
`rnspar.dat' is explained.  

A parameter file appears as follows

--
  120   24   48    8           : nrg, ntg, npg, nlg
   32   24   48    8           : nrf, ntf, npf, nlf
   30   40   JB   XZ           : nrf_deform, nrgin, NS_shape, EQ_point
   0.0d+00  1.25d+00   1.0d+04 : rgin, rgmid, rgout
#
 1800    y    n                : iter_max, sw_mass_iter, sw_art_deform
   0.4d-00   0.1d-00           : conv_gra, conv_ini
   0.4d-00   0.4d-00           : conv_den, conv_vep
   1D   3D   ccc               : indata_type,outdata_type,chrot,chgra,chope
   1.0e-06   1.0e-06           : eps, mass_eps
#
   20   -1                     : num_sol_seq, deform_par
  1.393927E-01  5.000000E-01   : emdc_ini, pinx
  8.396483E-02  7.559655E-02   : restmass_sph, gravmass_sph
  2.000000E-01  4.446856E-01   : MoverR_sph,   schwarz_radi_sph  (K=1 unit)
--


Most of the parameters are explained in the course notes.  
Functions of several flags are explained below.

NS_shape : flag to choose either an axi-symmetric or a tri-axial solution.  
           One of these symmetries is imposed on the matter variables.
NS_shape = JB  ->  a tri-axial configuration.      (JB taken from Jacobi) 
         = ML  ->  an axi-symmetric configuration. (ML taken from Maclaurin) 

(Note: Even if choosing JB, a solution may not necessarily become tri-axial. 
In a certain sub-space of parameters (M/R, axis ratio, EOS parameters, etc.), 
tri-axial solutions may not exist.  More often, there exist an axi-symmetric 
solution having the same parameters as a triaxial solution.  In this case, 
it depends on the initial guess of the iteration which solution is calculated.)

EQ_point : flag to choose the points for parameter equations.  
           Among three parameter eqs, two of them are at the center of the 
           star and the surface along the x-axis.  
EQ_point = XZ  ->  the third eq. is taken at the pole. 
         = XY  ->  the third eq. is taken at the surface along the y-axis
(Note: XZ may be good for any cases.  XY may be used to compute tri-axial
solution but iteration often failes so far.)

sw_mass_iter : iteration to calculate a solution with a given rest mass. 
sw_mass_iter = y  ->  iterate for a given rest mass. 
             = n  ->  do no iterate for a given rest mass. 

indata_type, outdata_type : data type for input and output data sets.
indata_type = 1D  ->  input data is read from `rns***_1D.ini'.
            = 3D  ->  input data is read from `rns***_3D.ini'.
outdata_type= 1D  ->  otuput data is written into `rns***_1D.las'.
            = 3D  ->  output data is written into `rns***_3D.las'.

num_sol_seq, deform_par : parameters to computate a sequence of solutions.
num_sol_seq = [integer#]  -> number of solutions to calculate.   

nrf_deform : the grid number where the third parameter equation is taken.
             "nrf_deform" is updated as "nrf_deform + deform_par".  
