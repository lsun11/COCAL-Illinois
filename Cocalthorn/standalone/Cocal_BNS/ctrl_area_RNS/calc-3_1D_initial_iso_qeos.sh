#! /usr/bin/tcsh -f
cp rnspar.dat  ../spherical_data/initial_isotropic_coord/work_area_1D_qeos/.
cp rnspar_surf.dat  ../spherical_data/initial_isotropic_coord/work_area_1D_qeos/.
cp qeos_parameter.dat   ../spherical_data/initial_isotropic_coord/work_area_1D_qeos/.
cp ../spherical_data/TOV_schwartzshild_coord/ovlas.dat  \
              ../spherical_data/initial_isotropic_coord/work_area_1D_qeos/.
cp ../spherical_data/TOV_schwartzshild_coord/ovphy_plot.dat  \
              ../spherical_data/initial_isotropic_coord/work_area_1D_qeos/.
cd ../spherical_data/initial_isotropic_coord/work_area_1D_qeos
./exe_Import_isotr_qeos
cd ../../../ctrl_area_RNS
