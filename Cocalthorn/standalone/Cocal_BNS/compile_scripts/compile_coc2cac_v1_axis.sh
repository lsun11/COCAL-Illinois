#! /bin/sh -f
cd ../code/Main_utility
gfortran -O3 -o exe_coc2cac coc2cac_v1_axis.f90

mv exe_coc2cac ../../../ID_BNS/.
rm *.mod
