#! /bin/sh -f
#pgf77 -Bstatic -fastsse -tp k8-64 -o 
cd ../code/Main_code
#/opt/intel/fce/9.1.040/bin/ifort -fast -o 
#gfortran -O3 -o exe_bh_test Main_poisson_bbh_test.f90
gfortran -O3 -o exe_BNS_CF  Main_corot_BNS_CF_mpt.f90

cp exe_BNS_CF ../../work_area_BNS/.
mv exe_BNS_CF ../../executable_files/.
rm *.mod
