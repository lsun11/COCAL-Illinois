#! /bin/sh -f

cd ../code/Main_utility/
gfortran -O3 -o exe_output_velocity output_fluid_velocity.f90
mv exe_output_velocity ../../executable_files/.
rm *.mod
