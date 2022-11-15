#!/bin/sh
#
grep -e 'rkijkij' -H */* -l   > file.list
grep -e 'rkijkij' -H */*/* -l >> file.list
#
LIST=./file.list
#
while read F1 
do
  sed -i 's/rkijkij/tfkijkij/g' ${F1}
  echo ${F1}
done < ${LIST}
#
rm ./file.list 
#
grep -e 'rkij' -H */* -l   > file.list
grep -e 'rkij' -H */*/* -l >> file.list
#
LIST=./file.list
#
while read F1 
do
  sed -i 's/rkij/tfkij/g' ${F1}
  echo ${F1}
done < ${LIST}
#
rm ./file.list 
#
grep -e 'rka' -H Subroutine/* -l   > file.list
LIST=./file.list
#
while read F1 
do
  sed -i 's/rka/tfka/g' ${F1}
  echo ${F1}
done < ${LIST}
#
rm ./file.list 
