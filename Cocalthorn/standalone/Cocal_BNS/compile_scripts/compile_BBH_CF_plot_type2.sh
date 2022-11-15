#! /bin/sh -f 
cd ../code/Main_utility 
gfortran -O3 -o exe_BBH_CF_plot interpolation_contour_potential_BBH_CF_type2.f90

cp exe_BBH_CF_plot ../../work_area_BBH/.
mv exe_BBH_CF_plot ../../executable_files/.
rm *.mod
