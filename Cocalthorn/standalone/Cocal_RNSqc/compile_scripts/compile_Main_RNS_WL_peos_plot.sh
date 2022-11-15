#! /bin/sh -f 
cd ../code/Main_utility 
gfortran -O3 -o exe_RNS_plot interpolation_contour_potential_RNS_WL.f90

cp exe_RNS_plot ../../work_area_RNS/.
mv exe_RNS_plot ../../executable_files/.
rm *.mod
