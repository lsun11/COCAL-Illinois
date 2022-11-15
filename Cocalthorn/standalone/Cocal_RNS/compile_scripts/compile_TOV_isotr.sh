#! /bin/sh -f
cd ../spherical_data/TOV_schwartzshild_coord/
#pgf77 -Bstatic -fastsse -tp k8-64 -o exe_TOV_peos TOV_peos.f90
#/opt/intel/fce/10.1.008/bin/ifort -fast -o exe_TOV_peos TOV_peos.f90
gfortran -O3 -o exe_TOV_peos TOV_peos_isotr.f90

cp exe_TOV_peos ../../executable_files/.
rm *.mod
