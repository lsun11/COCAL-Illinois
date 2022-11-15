#! /bin/sh -f
cd ../spherical_data/initial_isotropic_coord/code/
gfortran -O3 -o exe_Import_isotr_surf  Import_isotr_peos_surf.f90

cp exe_Import_isotr_surf  ../../../executable_files
cp exe_Import_isotr_surf  ../work_area_1D_peos
rm *.mod
