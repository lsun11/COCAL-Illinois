#! /usr/bin/tcsh -f

cd ../work_area_RNS

#mv rnsflu_1D.ini rnsflu_1D.ini.backup
#mv rnsgra_1D.ini rnsgra_1D.ini.backup
cp ../spherical_data/initial_isotropic_coord/work_area_1D_qeos/rnsflu_1D.las \
                                                               rnsflu_1D.ini
cp ../spherical_data/initial_isotropic_coord/work_area_1D_qeos/rnsgra_1D.las \
                                                               rnsgra_1D.ini
cd ../ctrl_area_RNS
