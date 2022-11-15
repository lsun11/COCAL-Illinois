#! /bin/sh -f 
cd ../code/Main_utility 
gfortran -O3 -o exe_MRNS_plot interpolation_contour_potential_MRNS.f90

cp exe_MRNS_plot ../../work_area_RNS/.
mv exe_MRNS_plot ../../executable_files/.
rm *.mod
