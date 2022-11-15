#!/bin/sh
#
rm tmp_include_*.list
rm tmp_file_*.list

grep -e 'use'  Main_code/Main_MagneticRNS_WL_peos.f90 > tmp_file_Module.list
grep -e 'call' Main_code/Main_MagneticRNS_WL_peos.f90 > tmp_file_Subroutine.list
#
LIST=./tmp_file_Module.list
#
while read F1 F2 F3
do
  echo "${F2}.f90" >> tmp_include_modulefiles_MRNS.list
done < ${LIST}
#
sed -i "s/,.f90/.f90/g" tmp_include_modulefiles_MRNS.list

LIST=./tmp_file_Subroutine.list
#
while read F1 F2 F3 F4 F5
ARRAY=( ${F1} ${F2} ${F3} ${F4} ${F5} )
do
  if test ${#ARRAY[*]} -eq 0 
  then 
    break
  else 
    for (( i = 0; i < ${#ARRAY[*]}; i++ ))
    {
      aaa=`expr match ${ARRAY[i]} call`
      if test ${aaa} -eq 4 
      then 
        echo "${ARRAY[i+1]}.f90" >> tmp_include_subroutines_MRNS.list
      fi
    }
  fi
done < ${LIST}
#
sed -i "s/(.*)//g" tmp_include_subroutines_MRNS.list
# rm ./file.list 

count=1
while [ $count -le 5 ]
do
echo "$count"
count=`expr $count + 1`
while read F1 F2 F3 F4 F5
LIST=./tmp_file_Subroutine.list
#


done

echo "end"
