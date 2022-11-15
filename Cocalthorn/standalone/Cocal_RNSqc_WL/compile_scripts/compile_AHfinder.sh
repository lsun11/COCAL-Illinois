#! /bin/sh -f
#pgf77 -Bstatic -fastsse -tp k8-64 -o 
filename_main=Main_AHfinder.f90
filename_exec=exe_AHfinder
cd ../code/Main_code
#/opt/intel/fce/9.1.040/bin/ifort -fast -o 
#gfortran -O3 -o exe_bh_test Main_poisson_bbh_test.f90
gfortran -O3 -o ${filename_exec} ${filename_main}

cp ${filename_exec} ../../work_area_BBH/.
mv ${filename_exec} ../../executable_files/.
rm *.mod
