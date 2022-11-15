#! /bin/sh -f 
cd ../code/Main_utility 
gfortran -O3 -o exe_BBH_CF_plot interpolation_contour_potential_BBH_CF.f90
mv exe_BBH_CF_plot ../../executable_files/.
