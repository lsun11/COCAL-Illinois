#! /bin/sh -f
#pgf77 -Bstatic -fastsse -tp k8-64 -o 
cd ../code/Main_code
#/opt/intel/fce/9.1.040/bin/ifort -fast -o 
gfortran -O3 -o exe_helm_test Main_helmholtz_bbh_test.f90

cp exe_helm_test ../../work_area_helmholtz_test/.
mv exe_helm_test ../../executable_files/.
rm *.mod
