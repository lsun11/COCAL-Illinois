#! /bin/sh -f

cd ../code/Main_utility/
gfortran -O3 -o exe_crust calc_crust_interface.f90
mv exe_crust ../../executable_files/.
rm *.mod
