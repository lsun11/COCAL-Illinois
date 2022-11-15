#! /usr/bin/tcsh -f
echo 'Data file is copied to INI_DATA directory  '

cd  ../spherical_data/initial_isotropic_coord/work_area_1D_peos

cp -i rnsflu_1D.las ../INI_data/INI_N05_C017.flu
cp -i rnsgra_1D.las ../INI_data/INI_N05_C017.gra
