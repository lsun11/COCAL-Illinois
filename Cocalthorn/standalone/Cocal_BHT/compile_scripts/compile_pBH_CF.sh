#! /bin/sh -f
#pgf77 -Bstatic -fastsse -tp k8-64 -o 
cd ../code/Main_code
#/opt/intel/fce/9.1.040/bin/ifort -fast -o 
#gfortran -O3 -o exe_bh_test Main_poisson_bbh_test.f90
gfortran -O3 -o exe_pBH_CF Main_pBH_CF.f90

cp exe_pBH_CF ../../work_area_BBH/.
mv exe_pBH_CF ../../executable_files/.
rm *.mod
