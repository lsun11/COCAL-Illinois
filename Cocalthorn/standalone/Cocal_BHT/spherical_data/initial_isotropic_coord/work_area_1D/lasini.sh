#! /usr/local/bin/tcsh -f
mv frgflu.ini frgflu.tmp
mv frggra.ini frggra.tmp
mv frgflu.las frgflu.ini
mv frggra.las frggra.ini
cp frgpar.dat frgflu.las
cp frgpar.dat frggra.las
