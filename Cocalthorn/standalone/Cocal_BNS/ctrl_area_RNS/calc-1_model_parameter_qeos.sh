#! /usr/bin/tcsh -f
cp rnspar.dat ../spherical_data/TOV_schwartzshild_coord/.
cp ovpar_qeos.dat ../spherical_data/TOV_schwartzshild_coord/.
cp qeos_parameter.dat ../spherical_data/TOV_schwartzshild_coord/.
cd ../spherical_data/TOV_schwartzshild_coord/
./exe_TOV_qeos
cd ../../ctrl_area_RNS
