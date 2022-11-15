#! /bin/sh -f
cd ../code/Main_utility
gfortran -O3 -o exe_plot interpolation_contour_potential_bbh_2pot_test_mpt.f90
mv exe_plot ../../executable_files/.
