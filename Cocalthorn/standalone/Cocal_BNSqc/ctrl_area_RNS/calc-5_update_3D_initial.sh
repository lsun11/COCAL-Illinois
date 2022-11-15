#! /usr/bin/tcsh -f

cd ../work_area_RNS
cp rnsflu_3D.ini rnsflu_3D.ini.tmp
cp rnsgra_3D.ini rnsgra_3D.ini.tmp
cp rnsgrids_3D.ini rnsgrids_3D.ini.tmp
cp rnsflu_3D.las rnsflu_3D.las.tmp
cp rnsgra_3D.las rnsgra_3D.las.tmp
cp rnsgrids_3D.las rnsgrids_3D.las.tmp
mv rnsflu_3D.las rnsflu_3D.ini
mv rnsgra_3D.las rnsgra_3D.ini
mv rnsgrids_3D.las rnsgrids_3D.ini
cp rnspar.dat rnsflu_3D.las
cp rnspar.dat rnsgra_3D.las 
cp rnspar.dat rnsgrids_3D.las
