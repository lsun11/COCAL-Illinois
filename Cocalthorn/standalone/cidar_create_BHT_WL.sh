# Last modified May 9, 2018
#
STARTTIME=$(date +%s)

#############################################################################
#  Select what kind of system to simulate:                                  #
#  BNS  : Binary neutron stars                                              #
#  BNSqc: Binary neutron stars with quark core                              #  
#  BBH  : Binary black holes                                                #
#  BHT  : Black hole torus system                                           #
#  RNS  : Rotating NS is conformal flat (CF) approximation                  #
#  RNSwl: Rotating NS is waveless formalism                                 #
#  RNSqc: Rotating NS in CF with quark core                                 #
#                                                                           #
compact_object=BHT
#                                                                           #
#############################################################################

original_code=Cocal_orig
new_dir_name=Cocal_${compact_object}
path_to_src=~/repositories/cocal_illinois/Cocalthorn
path_to_lib=~/repositories/cocal_illinois/lib
FCOMP=gfortran
#FCOMP=ifort
FFLAGS=
#-O2

if [ ! -d "${original_code}" ]; then
  echo "Directory ${original_code} does not exist."
  echo "${original_code} should have the new version of the standalone COCAL code...exiting"
  exit 1
fi

rm -rf ${new_dir_name}
cp -rfp ./Cocal_orig   ${new_dir_name}

echo "Preparing files..."
cd ${new_dir_name}/compile_scripts

if [ ${compact_object} == BNS ]; then
  echo "Compact object is BNS."
  sh select_weight_calc_midpoint_grav_th4th.sh
elif [ ${compact_object} == BNSqc ]; then                                        
  echo "Compact object is BNS with quark core."                                  
  sh select_weight_calc_midpoint_grav_th4th.sh                                   
  sh select_peos_qc.sh                                
elif [ ${compact_object} == BBH ]; then
  echo "Compact object is BBH."
  sh select_grid_r_bhex.sh
  sh select_grgrad_midpoint_r3rd_bhex.sh
  sh select_grgrad_gridpoint_4th_bhex.sh
  sh select_weight_calc_midpoint_grav_th4th.sh
elif [ ${compact_object} == BHT ]; then
  echo "Compact object is BHT."
  sh select_grgrad_midpoint_r3rd_bhex.sh
  sh select_grgrad_gridpoint_4th_bhex.sh
  sh select_weight_calc_midpoint_grav_th4th.sh
  sh select_set_read_23_15.sh
elif [ ${compact_object} == RNS ]; then
  echo "Compact object is RNS in conformal flat approximation."
  sh select_grgrad_midpoint_r3rd_NS.sh
  sh select_grgrad_gridpoint_4th_NS.sh
  sh select_weight_calc_midpoint_grav_th4th.sh
  sh select_reset_fluid_radius_le_1.sh
  sh select_correct_matter_source_out.sh
  sh select_calc_surface_quad.sh
elif [ ${compact_object} == RNSwl ]; then
  echo "Compact object is RNS in waveless formalism."
  sh select_grgrad_midpoint_r3rd_NS.sh
  sh select_grgrad_gridpoint_4th_NS.sh
  sh select_weight_calc_midpoint_grav_th4th.sh
  sh select_reset_fluid_radius_le_1.sh
  sh select_correct_matter_source_out.sh
  sh select_calc_surface_quad.sh
elif [ ${compact_object} == RNSqc ]; then
  echo "Compact object is RNS in CF formalism with quark core."
  sh select_grgrad_midpoint_r3rd_NS.sh
  sh select_grgrad_gridpoint_4th_NS.sh
  sh select_weight_calc_midpoint_grav_th4th.sh
  sh select_reset_fluid_radius_le_1.sh
  sh select_correct_matter_source_out.sh
  sh select_calc_surface_quad.sh
  sh select_peos_qc.sh
else 
  echo "Flag ${compact_object} should be one of {BNS,BNSqc,BBH,BHT,RNS,RNSwl,RNSqc}...exiting"
  exit 1
fi
sh cactus_eps_conflict.sh


echo "Copying files to" ${path_to_src}"/src/${new_dir_name}/"
cd ../code/Main_utility
/bin/rm -rf  ${path_to_src}/src/${new_dir_name}
/bin/mkdir   ${path_to_src}/src/${new_dir_name}
/bin/mkdir   ${path_to_src}/src/${new_dir_name}/Main_utility

/bin/cp -f  coc2pri_bht_WL.f90       ${path_to_src}/src/${new_dir_name}/Main_utility/

cd ..
/bin/cp -rf Analysis                 ${path_to_src}/src/${new_dir_name}/
/bin/cp -rf EOS                      ${path_to_src}/src/${new_dir_name}/
/bin/cp -rf Function                 ${path_to_src}/src/${new_dir_name}/
/bin/cp -rf Include_file             ${path_to_src}/src/${new_dir_name}/
/bin/cp -rf Module                   ${path_to_src}/src/${new_dir_name}/
/bin/cp -rf Module_interface         ${path_to_src}/src/${new_dir_name}/
/bin/cp -rf Module_mpatch            ${path_to_src}/src/${new_dir_name}/
/bin/cp -rf Subroutine               ${path_to_src}/src/${new_dir_name}/
/bin/cp -rf Subroutine_mpatch        ${path_to_src}/src/${new_dir_name}/

cd ${path_to_src}/src/${new_dir_name}/Main_utility/

echo "Compiling code..."
${FCOMP} ${FFLAGS} -c coc2pri_bht_WL.f90;   rm *.mod

echo "Creating library..."
ar rcvf libcocalbht.a coc2pri_bht_WL.o

rm coc2pri_bht_WL.o

cp libcocalbht.a ${path_to_lib}
echo "Copied library libcocalbht.a" to ${path_to_lib}
ENDTIME=$(date +%s)
echo "----------------------------------------------------------------"
echo "It took $(($ENDTIME - $STARTTIME)) seconds to complete this run."
