#! /bin/sh -f
cd ../code/Main_utility
gfortran -O3 -o exe_coc2cac   coc2cac_mrs_v1.f90

mv exe_coc2cac ../../../ID_TEST/.
rm *.mod