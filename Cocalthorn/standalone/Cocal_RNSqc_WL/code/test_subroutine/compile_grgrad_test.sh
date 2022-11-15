#! /bin/sh -f
#pgf77 -Bstatic -fastsse -tp k8-64 -o 
#gfortran -O3 -o exe_bh_test Main_poisson_bbh_test.f90
gfortran -O3 -o exe_test grgrad_4th_bhex_test.f90

rm *.mod
