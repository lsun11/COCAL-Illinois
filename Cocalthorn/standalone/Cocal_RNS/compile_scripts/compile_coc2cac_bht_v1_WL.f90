#! /bin/sh -f
cd ../code/Main_utility
gfortran -O3 -o exe_coc2cac  coc2cac_bht_v1_WL.f90

mv exe_coc2cac ../../../ID_TEST/.
rm *.mod
