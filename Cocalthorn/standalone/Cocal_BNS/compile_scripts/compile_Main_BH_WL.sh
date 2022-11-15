#! /bin/sh -f
#pgf77 -Bstatic -fastsse -tp k8-64 -o 
cd ../code/Main_code
#/opt/intel/fce/9.1.040/bin/ifort -fast -o 
gfortran -O3 -o exe_BH_WL Main_BH_WL.f90

cp exe_BH_WL ../../work_area_BH/.
mv exe_BH_WL ../../executable_files/.
rm *.mod
