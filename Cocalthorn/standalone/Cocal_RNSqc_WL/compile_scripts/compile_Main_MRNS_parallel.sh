#! /bin/sh -f
#pgf77 -Bstatic -fastsse -tp k8-64 -o 
cd ../code/Main_code
#/opt/intel/fce/9.1.040/bin/ifort -fast -o 
gfortran -O3 -fopenmp -o exe_MRNS_WL_peos Main_MagneticRNS_WL_peos.f90

cp exe_MRNS_WL_peos ../../work_area_RNS/.
mv exe_MRNS_WL_peos ../../executable_files/.
rm *.mod
