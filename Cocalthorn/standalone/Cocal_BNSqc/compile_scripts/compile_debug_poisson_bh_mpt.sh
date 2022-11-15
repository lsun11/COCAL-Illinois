#! /bin/sh -f
#pgf77 -Bstatic -fastsse -tp k8-64 -o 
cd ../code/Main_code
#/opt/intel/fce/9.1.040/bin/ifort -fast -o 
gfortran -O3 -o exe_mpt_test Main_poisson_bbh_test_mpt.f90

cp exe_mpt_test ../../work_area_poisson_test/.
mv exe_mpt_test ../../executable_files/.
rm *.mod
