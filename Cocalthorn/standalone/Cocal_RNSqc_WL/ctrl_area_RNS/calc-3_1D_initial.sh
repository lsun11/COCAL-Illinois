#! /usr/bin/tcsh -f
echo 'CHOOSE INITIAL DATA CORRECTLY '
cp rnspar.dat ../spherical_data/initial_isotropic_coord/work_area_1D_peos/.
cp peos_parameter.dat \
              ../spherical_data/initial_isotropic_coord/work_area_1D_peos/.
cd ../spherical_data/initial_isotropic_coord/work_area_1D_peos
./exe_Main_1Dini_peos
cd ../../../ctrl_area_RNS
