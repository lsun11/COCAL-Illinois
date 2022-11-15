#! /usr/bin/tcsh -f
echo 'Once again'

cd  ../spherical_data/initial_isotropic_coord/work_area_1D_peos
#sh nxtini.sh
cp rnsflu_1D.ini ../../../work_area/rnsflu_1D.tmp
cp rnsgra_1D.ini ../../../work_area/rnsgra_1D.tmp

cp rnsflu_1D.nxt ../../../work_area/rnsflu_1D.ini
cp rnsgra_1D.nxt ../../../work_area/rnsgra_1D.ini
./exe_Main_1Dini_peos

