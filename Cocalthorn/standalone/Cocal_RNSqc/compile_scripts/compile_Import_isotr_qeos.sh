#! /bin/sh -f
cd ../spherical_data/initial_isotropic_coord/code/
gfortran -O3 -o exe_Import_isotr_qeos  Import_isotr_qeos.f90

cp exe_Import_isotr_qeos  ../../../executable_files
cp exe_Import_isotr_qeos  ../work_area_1D_qeos
rm *.mod
