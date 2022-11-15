#! /bin/sh -f
#pgf77 -Bstatic -fastsse -tp k8-64 -o 
cd ../code/Main_utility
#/opt/intel/fce/9.1.040/bin/ifort -fast -o 
gfortran -O3 -o exe_BHT_WL_plot  interpolation_contour_potential_BHT_WL.f90

cp exe_BHT_WL_plot ../../work_area_BHT/.
mv exe_BHT_WL_plot ../../executable_files/.
rm *.mod
