#! /usr/bin/tcsh -f
#pgf77 -Bstatic -fastsse -tp k8-64 -o 
cd ../code/Main_code
#pgf77 -Bstatic -fastsse -tp k8-64 -o exe_Main_1Dini_peos Main_1Dini_peos.f90
#/opt/intel/fce/10.1.008/bin/ifort -fast -o exe_Main_1Dini_peos Main_1Dini_peos.f90
gfortran -O3 -o exe_mpt_test Main_multipatch_test.f90

cp exe_mpt_test ../../executable_files
cp exe_mpt_test ../../work_area_poisson_test
rm *.mod
