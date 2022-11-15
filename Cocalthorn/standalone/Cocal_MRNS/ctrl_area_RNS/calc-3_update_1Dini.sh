#! /usr/bin/tcsh -f
echo 'Update 1D initial for 1D code'

cd  ../spherical_data/initial_isotropic_coord/work_area_1D_peos

mv rnsflu_1D.ini rnsflu_1D.tmp
mv rnsgra_1D.ini rnsgra_1D.tmp
cp rnsflu_1D.las rnsflu_1D.ini
cp rnsgra_1D.las rnsgra_1D.ini
