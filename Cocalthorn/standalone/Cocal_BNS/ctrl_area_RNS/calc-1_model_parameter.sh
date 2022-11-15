#! /usr/bin/tcsh -f
cp rnspar.dat ../spherical_data/TOV_schwartzshild_coord/.
cp ovpar_peos.dat ../spherical_data/TOV_schwartzshild_coord/.
cp peos_parameter.dat ../spherical_data/TOV_schwartzshild_coord/.
cd ../spherical_data/TOV_schwartzshild_coord/
./exe_TOV_peos
cd ../../ctrl_area_RNS
