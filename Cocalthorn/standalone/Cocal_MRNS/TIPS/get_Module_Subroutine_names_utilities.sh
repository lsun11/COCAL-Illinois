#!/bin/sh
#
# TAG=MRNS
# filename_Maincode=./Main_code/Main_MagneticRNS_WL_peos.f90
### TAG=AHfinder
### filename_Maincode=./Main_code/Main_AHfinder.f90
# TAG=BBH_CF
# filename_Maincode=./Main_code/Main_BBH_CF.f90
# TAG=RNS_CF_peos
# filename_Maincode=./Main_code/Main_RNS_CF_peos.f90
#Plot util RNS_WL TAG=RNS_WL_peos_plot
#Plot util RNS_WL filename_Maincode=./Main_utility/interpolation_contour_potential_RNS_WL.f90
#Plot util # TAG=RNS_CF_peos_plot
#Plot util # filename_Maincode=./Main_utility/interpolation_main.f90
#
#Plot util # 
TAG=plot_MRNS
#Plot util # 
filename_Maincode=./Main_utility/interpolation_contour_potential_MRNS.f90
#
#  TAG=plot_BBH_CF_3mpt
#  filename_Maincode=./Main_utility/interpolation_contour_potential_BBH_CF_3mpt.f90
#
 filename_module=./tmp_include_modulefiles_${TAG}.list_
 filename_subrou=./tmp_include_subroutines_${TAG}.list_
# --
dir_name_module=./Module/
dir_name_subrou=./Subroutine/
dir_name_module_ana=./Analysis/Module/
dir_name_subrou_ana=./Analysis/Subroutine/
###
incfile_module=./include_modulefiles_${TAG}.f90
incfile_interf=./include_interface_modulefiles_${TAG}.f90
incfile_subrou=./include_subroutines_${TAG}.f90
incfile_module_ana=./include_modulefiles_analysis_${TAG}.f90
incfile_interf_ana=./include_interface_modulefiles_analysis_${TAG}.f90
incfile_subrou_ana=./include_subroutines_analysis_${TAG}.f90
###
incdir_module=../Module/
incdir_interf=../Module_interface/
incdir_subrou=../Subroutine/
incdir_module_ana=../Analysis/Module/
incdir_interf_ana=../Analysis/Module/
incdir_subrou_ana=../Analysis/Subroutine/
# --
#Plot util #dir_name_module=./Analysis/Module/
#Plot util #dir_name_subrou=./Analysis/Subroutine/
#Plot util #incfile_module=./include_modulefiles_analysis_${TAG}.f90
#Plot util #incfile_interf=./include_interface_modulefiles_analysis_${TAG}.f90
#Plot util #incfile_subrou=./include_subroutines_analysis_${TAG}.f90
#Plot util #incdir_module=../Analysis/Module/
#Plot util #incdir_interf=../Module_interface/
#Plot util #incdir_subrou=../Analysis/Subroutine/
# --
#   
rm tmp_*.list*
rm tmp_*.list_*
rm tmp_mod.list
rm tmp_sub.list

name_clean_module(){
# Read input file name, output file name
#
# For modules
#
LIST_mod=./$1
while read F1 F2 F3
ARRAY=( ${F1} ${F2} ${F3} )
do
  if test ${#ARRAY[*]} -eq 0 
  then 
    break
  else 
    for (( i = 0; i < ${#ARRAY[*]}; i++ ))
    {
      if test "${ARRAY[i]}" = "use" 
      then 
        echo "${ARRAY[i+1]}.f90" >> $2
      fi
    }
  fi
done < ${LIST_mod}
#
sed -i "s/,.f90/.f90/g" $2
}
#
name_clean_subrou(){
# Read input file name, output file name
#
# For subroutines
#
LIST_sub=./$1
while read F1 F2 F3 F4 F5 F6 F7
ARRAY=( ${F1} ${F2} ${F3} ${F4} ${F5} ${F6} ${F7} )
do
  if test ${#ARRAY[*]} -eq 0 
  then 
    break
  else 
    for (( i = 0; i < ${#ARRAY[*]}; i++ ))
    {
      if test "${ARRAY[i]}" = "call" 
      then 
        echo "${ARRAY[i+1]}" >> $2
      fi
    }
  fi
done < ${LIST_sub}
#
sed -i "s/(.*//g" $2
sed -i "s/$/\.f90/g" $2
}
#
get_name_module_subroutine(){
# get lines contains use and call
# In order of input file name, 
#            output file name for module, 
#            output file name for subroutine
#
  grep -e 'use'  $1 >> $2
  grep -e 'call' $1 >> $3
}
sort_and_remove_redundancy(){
  cp $1 $1tmp ; rm $1 
  sort $1tmp | uniq > $1
  rm $1tmp
}
#

count=1
nest=7
while [ ${count} -le ${nest} ]
if test ${count} -eq `expr ${nest} + 1`
then 
  break
fi
do
  echo "$count"
  if test ${count} -eq 1
  then
    A1=${filename_Maincode} 
    A2=${filename_module}${count}
    A3=${filename_subrou}${count}
    get_name_module_subroutine ${A1} tmp_mod.list tmp_sub.list
  else
    count_m1=`expr $count - 1`
    N1=${filename_module}${count_m1}
    LIST=${N1}
    while read G1
    do 
      A1=${dir_name_module}${G1}
      A2=${filename_module}${count}
      A3=${filename_subrou}${count}
      get_name_module_subroutine ${A1} tmp_mod.list tmp_sub.list
    done < ${LIST}
    N1=${filename_module}${count_m1}
    LIST=${N1}
    while read G1
    do 
      A1=${dir_name_module_ana}${G1}
      A2=${filename_module}${count}
      A3=${filename_subrou}${count}
      get_name_module_subroutine ${A1} tmp_mod.list tmp_sub.list
    done < ${LIST}
    N1=${filename_subrou}${count_m1}
    LIST=${N1}
    while read G1
    do 
      A1=${dir_name_subrou}${G1}
      A2=${filename_module}${count}
      A3=${filename_subrou}${count}
      get_name_module_subroutine ${A1} tmp_mod.list tmp_sub.list
    done < ${LIST}
    N1=${filename_subrou}${count_m1}
    LIST=${N1}
    while read G1
    do 
      A1=${dir_name_subrou_ana}${G1}
      A2=${filename_module}${count}
      A3=${filename_subrou}${count}
      get_name_module_subroutine ${A1} tmp_mod.list tmp_sub.list
    done < ${LIST}
  fi
  echo "Cleaning"
  name_clean_module tmp_mod.list ${A2}
  name_clean_subrou tmp_sub.list ${A3}
  rm tmp_mod.list; rm tmp_sub.list
  echo "Sorting"
  sort_and_remove_redundancy ${A2}
  sort_and_remove_redundancy ${A3}
  count=`expr $count + 1`
done
#
# Merge files
# 
rm ${incfile_module}
rm ${incfile_interf}
rm ${incfile_subrou}
rm ${incfile_module_ana}
rm ${incfile_subrou_ana}
#
count=1
while [ ${count} -le ${nest} ]
if test ${count} -eq `expr ${nest} + 1`
then 
  break
fi
do
  A2=${filename_module}${count}
  A3=${filename_subrou}${count}
  cat ${A2} >> ${incfile_module}
  cat ${A3} >> ${incfile_subrou}
  rm ${A2} ; rm ${A3}
  count=`expr $count + 1`
done
#Remove redundancy
echo "Remove redundancy"
sort_and_remove_redundancy ${incfile_module}
sort_and_remove_redundancy ${incfile_subrou}
grep -e "interface_"   ${incfile_module} > ${incfile_interf}
sed -i "/interface_/d" ${incfile_module}
sed -i "/flusw\.f90/d"    ${incfile_module}
#sed -i "/use\.f90/d"    ${incfile_module}
#sed -i "s!Module/\.f90!delete_this_line!g"    ${incfile_module}
#sed -i "/delete_this_line/d"    ${incfile_module}
sed -i "/alloc_array/d"  ${incfile_subrou}
sed -i "/copy_array/d"   ${incfile_subrou}

sed -i "/alloc_int/d"  ${incfile_subrou}
sed -i "/allocate_hgfn/d"  ${incfile_subrou}
sed -i "/allocate_legendre/d"  ${incfile_subrou}
sed -i "/allocate_vector_/d"  ${incfile_subrou}
sed -i "/allocate_weight_/d"  ${incfile_subrou}
sed -i "/calc_parameter_binary_excision/d"        ${incfile_subrou}
sed -i "/calc_hgfn/d"        ${incfile_subrou}
sed -i "/^grid_/d"        ${incfile_subrou}
sed -i "/^peos_/d"        ${incfile_subrou}
#cp ${incfile_subrou} baka_0
sed -i "/trig_grav_/d"   ${incfile_subrou}
#cp ${incfile_subrou} baka_1
sed -i "/weight_calc_/d" ${incfile_subrou}
#cp ${incfile_subrou} baka_2
sed -i "/calc_weight_/d" ${incfile_subrou}
#cp ${incfile_subrou} baka_3
sed -i "/choose_formulation/d"   ${incfile_subrou}
sed -i "/hegr4/d"   ${incfile_subrou}
sed -i "/hrhafn/d"   ${incfile_subrou}
sed -i "/legendre\.f90/d"   ${incfile_subrou}
sed -i "/legendre_theta\.f90/d"   ${incfile_subrou}
sed -i "/read_parameter\.f90/d"   ${incfile_subrou}
sed -i "/read_parameter_binary_excision\.f90/d"   ${incfile_subrou}
sed -i "/read_parameter_cartesian\.f90/d"   ${incfile_subrou}
sed -i "/subroutines\./d" ${incfile_subrou}

ls ${dir_name_module_ana} > tmp_module_ana
ls ${dir_name_subrou_ana} > tmp_subrou_ana
#cp  ${incfile_module} test1
#cp  ${incfile_interf} test2
#cp  ${incfile_subrou} test3
comm -12 tmp_module_ana  ${incfile_module} > ${incfile_module_ana}
comm -12 tmp_module_ana  ${incfile_interf} > ${incfile_interf_ana}
comm -12 tmp_subrou_ana  ${incfile_subrou} > ${incfile_subrou_ana}
comm -13 tmp_module_ana  ${incfile_module} > tmp_module
comm -13 tmp_module_ana  ${incfile_interf} > tmp_interf
comm -13 tmp_subrou_ana  ${incfile_subrou} > tmp_subrou
cp tmp_module ${incfile_module} 
cp tmp_interf ${incfile_interf} 
cp tmp_subrou ${incfile_subrou} 
rm tmp_module*; rm tmp_interf*; rm tmp_subrou*

#Add directory names
echo "Add directory names"
sed -i "s!^!include '${incdir_module}!" ${incfile_module}
sed -i "s!^!include '${incdir_interf}!" ${incfile_interf}
sed -i "s!^!include '${incdir_subrou}!" ${incfile_subrou}
sed -i "s!^!include '${incdir_module_ana}!" ${incfile_module_ana}
sed -i "s!^!include '${incdir_interf_ana}!" ${incfile_interf_ana}
sed -i "s!^!include '${incdir_subrou_ana}!" ${incfile_subrou_ana}
sed -i "s/\$/'/" ${incfile_module}
sed -i "s/\$/'/" ${incfile_interf}
sed -i "s/\$/'/" ${incfile_subrou}
sed -i "s/\$/'/" ${incfile_module_ana}
sed -i "s/\$/'/" ${incfile_interf_ana}
sed -i "s/\$/'/" ${incfile_subrou_ana}

echo "end"
