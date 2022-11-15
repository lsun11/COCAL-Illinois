#! /bin/sh -f
#pgf77 -Bstatic -fastsse -tp k8-64 -o 
cd ../code/Main_code
#/opt/intel/fce/9.1.040/bin/ifort -fast -o 
#gfortran -O3 -o exe_bh_test Main_poisson_bbh_test.f90
#gfortran -std=gnu -O3 -o exe_BBH_CF Main_BBH_CF_circ.f90
gfortran -O3 -o exe_BBH_CF Main_BBH_CF_circ_spin.f90

cp exe_BBH_CF ../../work_area_BBH/.
mv exe_BBH_CF ../../executable_files/.
rm *.mod
