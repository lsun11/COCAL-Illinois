STARTTIME=$(date +%s)

######################## Create the TEST ID directory ############################
#
ID_dir_name="/home/antonios/Data/MY/MYPYTHON/MRNS/work_area_RNS_M1.5_SD3"
#ID_dir_name="/u/sciteam/tsokaros/COCAL_ID/spHg2.5_SLy_231_2cor/work_area_BNS"
#ID_dir_name="/home/astro/shared/Cocal_ID/K123_ir_15vs15_Hs2d_45km"
#ID_dir_name="/home/bruno/scratch/initial_data/Cocal_ID/K123_ir_15vs15_Hs2d_45km"

sed -i "s/d/e/g" ${ID_dir_name}/peos_parameter*.dat
sed -i "s/E/e/g" ${ID_dir_name}/peos_parameter*.dat

if [ ! -d "ID_TEST" ]; then
  echo "Create folder ID_TEST..."
  mkdir ID_TEST
else
  echo "Remove existing folder ID_TEST and create a new one..."
  rm -rf  ID_TEST
  mkdir ID_TEST
#  exit 1
fi
#
##################################################################################

new_dir_name=tmpdir_eraseme

rm -rf ${new_dir_name}
cp -rfp ./Cocal_orig   ${new_dir_name}
cd ${new_dir_name}

cd compile_scripts
##################################################################################
#sh compile_coc2cac.sh             # Original coc2cac. Only for irrotational.
#sh compile_coc2cac_co.sh          # Only for corotating.
#sh compile_coc2cac_ir.sh          # Only for irrotational.
#sh compile_coc2cac_sp.sh          # Only for spinning.
##################################################################################

#For binary neutron stars
#echo "compile_coc2cac_ini.sh..."
#sh compile_coc2cac_ini.sh      

#For binary neutron stars with quark core
#echo "compile_coc2cac_ini.sh..."
#sh select_peos_qc.sh
#sh compile_coc2cac_ini.sh      

#For rotating star in CF formalism
#echo "compile_coc2cac_rs_v1.sh..."
#sh compile_coc2cac_rs_v1.sh  

#For rotating star in WL formalism
#echo "compile_coc2cac_rs_v1_WL.sh..."
#sh compile_coc2cac_rs_v1_WL.sh   

#For magnetized rotating star
echo "compile_coc2cac_mrs_v1.sh..."
sh compile_coc2cac_mrs_v1.sh

#For BHT in WL formalism
#echo "compile_coc2cac_bht_v1_WL.sh..."
#sh select_set_read_23_15.sh       # By default is reading 20.12 for RNS WL
#sh compile_coc2cac_bht_v1_WL.f90  


cd ../../ID_TEST
for i in `/bin/ls $ID_dir_name `; do  ln -s ${ID_dir_name}/$i ;  done
rm -rf *.plt  *.f90  exe_BNS*  exe_TOV*  exe_Import*  exe_BHT*

./exe_coc2cac
#coc2cac_pid=`echo $!`
#echo $! > coc2cac_pid.txt

cd ..  # Now we are in Coc2Cac
cd ..

ENDTIME=$(date +%s)
echo "----------------------------------------------------------------"
echo "It took $(($ENDTIME - $STARTTIME)) seconds to complete this run."

