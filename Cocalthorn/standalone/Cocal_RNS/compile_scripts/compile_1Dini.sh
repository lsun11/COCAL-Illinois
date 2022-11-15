#! /bin/sh -f
#pgf77 -Bstatic -fastsse -tp k8-64 -o 
cd ../spherical_data/initial_isotropic_coord/code/
#pgf77 -Bstatic -fastsse -tp k8-64 -o exe_Main_1Dini_peos Main_1Dini_peos.f90
#/opt/intel/fce/10.1.008/bin/ifort -fast -o exe_Main_1Dini_peos Main_1Dini_peos.f90
gfortran -O3 -o exe_Main_1Dini_peos Main_1Dini_peos.f90

cp exe_Main_1Dini_peos ../../../executable_files
cp exe_Main_1Dini_peos ../work_area_1D_peos
rm *.mod
