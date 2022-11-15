#! /bin/sh -f 
cd ../code/Main_utility 
gfortran -O3 -o exe_RQS_plot interpolation_main_qeos.f90

cp exe_RQS_plot ../../work_area_RNS/.
mv exe_RQS_plot ../../executable_files/.
rm *.mod
