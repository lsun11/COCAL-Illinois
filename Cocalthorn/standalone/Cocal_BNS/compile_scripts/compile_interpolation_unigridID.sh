#! /bin/sh -f

cd ../code/Main_utility/
gfortran -O3 -o exe_interpolation_unigridID interpolation_unigridID.f90
mv exe_interpolation_unigridID ../../executable_files/.
rm *.mod
