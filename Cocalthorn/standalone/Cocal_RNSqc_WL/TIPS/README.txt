Cocal - Compact Object Calculator project

Cocal is a minimal code set for computing 
relativistic compact objects in (quasi-)equilibrium, 
including rotating neutron stars, 
rotating neutron stars with magnetic fields, 
binary neutron stars, 
binary neutron stars with magnetic fields, 
black hole - neutron star binaries, 
black hole - neutron star binaries with magnetic fields, 
binary black holes. 
It is aimed to provide a minimal code set for 
computing numerical solutions of these objects, 
in which used are 2nd order (or higher) finite 
difference scheme on spherical coordinates, 
numerical quadratures of multipole expansions 
of Green's function, and self consistent iterations.  
Cocal is written in fortran90 language.  


----

Basic structure of directories is as follows.

/Cocal -- /code ==> Main_code, Modules, Subroutines, Functions, EOS, etc.
       |
       -- /compile_scripts ==> shell scripts to compile codes
       |
       -- /executable_files ==> a storage for executable files
       |
       -- /parameters_files ==> a storage for parameter files
       |
       -- /ctrl_area ==> shell scripts, parameter files.
       |
       -- /work_area ==> input and output data, physical quantities, etc.
       |
       -- /spherical_data --(*) ==> code and data for computing spherical stars


(*) -- /TOV_schwartzshild_coord ==> TOV code, data 
    |
    -- /initial_isotropic_coord --- /code ==> 1D code on isotropic coordinate
                                 |
                                 -- /work_area_1D ==> work area for 1D code
                                 |
                                 -- /INI_data ==> data files for 1D solutions


--

Assuming the code being fully debugged, whole calculations can be 
operated using the shell scripts in the directory `/Cocal/ctrl_area'.
It is recommended to copy this directory, modify contents of 
parameter files, and then run the scripts

As a matter of fact, whenever those scripts are issued, 
it is recommended to check their contents on the console beforehand.  

The main code is executed in the directory `/Cocal/work_area'.  
After a successful run, all data for the solutions will be over 
written in the files in this directory.  

The codes may be compiled and run following the steps below; 

0. Compile all codes.  

   (1) Modify the compiler command of the shell scripts in 
       `/Cocal/compile_scripts' according to the fortran90 compiler 
       on your system, then issue those scripts.  To compile 
       frequently used codes, run a script 'update_all.sh' in this 
       directory.  

The followings are done in `/Cocal/ctrl_area' directory.

1. Set the grid and EOS parameters.

   (1) Modify the grid parameters in `rnspar.dat' if it is necessary.
       One can copy typical parameters (e.g. the parameters for the 
       resolutions) from a sample file such as, `rnspar.dat.type1-large'. 
       Sample files are found in `/Cocal/parameters/rns_code_parameters/'.

   (2) Set the EOS parameters in `peos_parameters.dat'.  Sample EOS 
       parameters are found in `/Cocal/parameters/peos_parameters/'.

2. Compute the model parameters.  

   In the output of the Cocal code, certain quantities are normalized 
   by the gravitational mass (or other quantities) of the spherical 
   solution having the same rest mass.  Using the TOV solver, such 
   quantities of corresponding spherical solutions are calculated.  

   To do this, follow the steps:

   (1) Choose a value to fix which may be the compactness M/R, 
       the rest mass, or the gravitational mass, which is controled by 
       the parameters `itype' and `iter_mode' in `ovpar_peos.dat' file.  

       Edit those numbers in the file `rnspar.dat', either 
       one of  `MoverR_sph', `restmass_sph', or `gravmass_sph'.  
       Set an approximate value of rho at the center in cgs unit
       in 'peos_parameter.dat' file, which is `rhoini_cgs'.  

       Note 1. rhoini_cgs is translated to p/rho 'emdc_ini' in 
               G=c=Msol=1 unit.  `pinx' in `rnspar.dat' file is 
               not used in parametrized EOS (which stands a polytopic 
               index of the polytrope code).  

       Note 2. that the fixed formats are used in the parameter files.  

       Note 3. Other parameters in `rnspar.dat' are nothing to do 
               with TOV code.  Grid spacing of this TOV solver is 
               taken from the data file in `ovpar_peos.dat'.  
               It is not necessary to change the resolution of TOV 
               solver unless one wishes to have the higher precision.


   (2) Run a shell script 

       `calc-1_model_parameter.sh'.  

       This script copies `rnspar.dat' and `peos_parameter.dat' to 
       the TOV directory, and run the TOV solver.  The result may be 
       used as model parameters for RNS calculations in an output file 
       `rnspar_add.dat', which has the same format as the last three 
       lines of `rnspar.dat'.  Physical quantities of TOV solution 
       are written in `ovphy.dat'.  

2. Run a shell script 

   `calc-2_update_rnspar.sh'

    This script replaces the last 3 lines of `rnspar.dat', calculated 
    in the step 1 above.  


3. Compute a spherical solution used for the initial guess for rotating 
   star computations, using

  `Main_1Dini_peos.f90'

   This code calculate a spherical solution in isotropic coordinate, 
   whose parameter is specified in `rnspar.dat'. The solution is used 
   for the inital guess of the rotating star code.

   (1) If the initial guess for 1D code is largely different from 
       the parameter of your trial, choose an example for the initial 
       data from `INI_data' directory.  
       One can MODIFY the shell script 'calc-3_choose_1Ddata.sh'
       and run it to update the initial guess for the 1D code.  

   (2) Run a shell script

       `calc-3_1D_initial.sh'. 

       With this shell script, `rnspar.dat' and `peos_parameter.dat' 
       files are copied to the `/work_area_1D' directory, and the 1D 
       code is executed.  

       (2-1) If the iteration of the 1D code gives "NOT A NUMBER", 
             modify the parameters and/or the initial guess so that 
	     the difference between parameters of the initial guess 
	     and those of new calculation to be smaller.  

       (2-2) If the iteration of (2) didn't converge, run a script

             `calc-3_repeat_1D_initial.sh'

              to repeat the iteration.  

       (2-3) When the iteration converged, one can MODIFY the shell 
	     script 'calc-4_update_1D_initial.sh' to keep the data
             in `Cocal/spherical_data/initial_isotropic_coord/INI_data', 
	     if necessary.  

4. Run a shell script 

   `calc-4_update_1D_initial.sh'.

    This script replace the inital data for the rotating star code 
    by the new data.  If you already have the data of a spherical solution 
    in `Cocal/spherical_data/initial_isotropic_coord/INI_data' directory, 
    one can skip the step 3 above, and MODIFY `calc-4_update_1D_initial.sh' 
    just to copy the data from this directory.  

5. Run the rotating star code.

   `calc-5_compute_RNS.sh'

   Once the iteration is converged, the solution can be used 
   as the initial guess for the next computation.  To copy the 
   data files of converged data to that of the initial data, 
   one can use a shell script `calc-5_update_3D_initial.sh'.  

(Note: Since these are not finalized, the code and structure of 
directories may be modified in the future.  The other code woudn't 
be modified unless we modify the definition of grid parameters. ) 


--
Computing a sequence of solutions.

In the current version of the code, when a sequece of solutions 
for compact stars is calculated, the rest mass of each solution 
is set to be constant.  However, for a better presentation of 
the results, we specify the sequence by the compactness M/R of 
a spherical solution having the same rest mass.  
We first choose a certain value of the compactness, and then 
calculate the corresponding rest mass of a spherical solution.  
This value of the rest mass is used to construct a constant 
rest mass sequnece of rotating stars.  

Explained below, there are 3 codes for the doing these calculations 
which are placed in different directories.  

Those codes are :

1. `spherical_data/TOV_schwartzshild_coord/TOV_peos.f90'
    computes a TOV solution, which gives the rest mass and M/R.  

2. `spherical_data/initial_isotropic_coord/Main_1Dini_peos.f90'
    computes a spherical solution in the isotropic coordinate, 
    which is used for the initial data for the deformed stars.

3. `code/Main_code/Main_RNS_CF_peos.f90'
    computes a rotating neutron star.  

