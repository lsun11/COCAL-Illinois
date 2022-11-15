#! /usr/bin/tcsh -f
mv rnspar.dat rnspar.dat.backup
head -12 rnspar.dat.backup > rnspar.dat
cat ../spherical_data/TOV_schwartzshild_coord/rnspar_add.dat >> rnspar.dat
