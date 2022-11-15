#! /usr/bin/tcsh -f
cp rnspar.dat  ../spherical_data/initial_isotropic_coord/work_area_1D_peos/.
cp rnspar_surf.dat  ../spherical_data/initial_isotropic_coord/work_area_1D_peos/.
cp peos_parameter.dat   ../spherical_data/initial_isotropic_coord/work_area_1D_peos/.
cp ../spherical_data/TOV_schwartzshild_coord/ovlas.dat  \
              ../spherical_data/initial_isotropic_coord/work_area_1D_peos/.
cp ../spherical_data/TOV_schwartzshild_coord/ovphy_plot.dat  \
              ../spherical_data/initial_isotropic_coord/work_area_1D_peos/.
cd ../spherical_data/initial_isotropic_coord/work_area_1D_peos
./exe_Import_isotr_surf
cd ../../../ctrl_area_RNS
