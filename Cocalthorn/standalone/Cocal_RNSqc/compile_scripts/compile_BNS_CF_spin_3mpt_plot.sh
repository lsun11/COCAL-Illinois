#! /bin/sh -f
cd ../code/Main_utility
gfortran -O3 -o exe_BNS_CF_plot interpolation_contour_potential_spin_BNS_CF_3mpt.f90

cp exe_BNS_CF_plot ../../work_area_BNS/.
mv exe_BNS_CF_plot ../../executable_files/.
rm *.mod
