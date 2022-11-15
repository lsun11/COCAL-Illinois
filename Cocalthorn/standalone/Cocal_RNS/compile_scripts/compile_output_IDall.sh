#! /bin/sh -f

cd ../code/Main_utility/
gfortran -O3 -o exe_output_IDall output_IDall.f90
mv exe_output_IDall ../../executable_files/.
rm *.mod
