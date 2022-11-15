#! /bin/sh -f
#pgf77 -Bstatic -fastsse -tp k8-64 -o 
cd ../code/Main_code
#/opt/intel/fce/9.1.040/bin/ifort -fast -o 
gfortran -O3 -o exe_RNS_CF_peos Main_RNS_CF_peos.f90

cp exe_RNS_CF_peos ../../work_area_RNS/.
mv exe_RNS_CF_peos ../../executable_files/.
rm *.mod
