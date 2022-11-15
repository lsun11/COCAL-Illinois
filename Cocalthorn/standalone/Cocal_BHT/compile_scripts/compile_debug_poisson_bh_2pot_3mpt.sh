#! /bin/sh -f
#pgf77 -Bstatic -fastsse -tp k8-64 -o 
cd ../code/Main_code
#/opt/intel/fce/9.1.040/bin/ifort -fast -o 
#gfortran -O3 -o exe_bh_test Main_poisson_bbh_test.f90
gfortran -O3 -o exe_bh_test Main_poisson_bbh_2pot_test_3mpt.f90

cp exe_bh_test ../../work_area_poisson_test/.
mv exe_bh_test ../../executable_files/.
rm *.mod
