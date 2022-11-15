#! /bin/sh -f
#pgf77 -Bstatic -fastsse -tp k8-64 -o 
cd ../code/Main_code
#/opt/intel/fce/9.1.040/bin/ifort -fast -o 
gfortran -O3 -o exe_RQS_CF_qeos Main_RQS_CF_qeos.f90

cp exe_RQS_CF_qeos ../../work_area_RNS/.
mv exe_RQS_CF_qeos ../../executable_files/.
rm *.mod
