#! /bin/sh -f
cd ../spherical_data/initial_isotropic_coord/code/
gfortran -O3 -o exe_Import_isotr  Import_isotr_peos.f90

cp exe_Import_isotr  ../../../executable_files
cp exe_Import_isotr  ../work_area_1D_peos
rm *.mod
