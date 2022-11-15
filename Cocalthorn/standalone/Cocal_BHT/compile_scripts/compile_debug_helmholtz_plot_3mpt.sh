#! /bin/sh -f
cd ../code/Main_utility
#gfortran -O3 -o exe_plot interpolation_contour_potential_bbh_test.f90
#gfortran -O3 -o exe_plot interpolation_contour_wave_test.f90
gfortran -O3 -o exe_helm_plot interpolation_contour_wave_binary_test_mpt.f90
mv exe_helm_plot ../../executable_files/.
