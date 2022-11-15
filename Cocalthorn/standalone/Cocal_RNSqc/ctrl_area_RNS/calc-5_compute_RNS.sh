#! /usr/bin/tcsh -f
cp rnspar.dat ../work_area_RNS/.
cp peos_parameter.dat ../work_area_RNS/.
cd ../work_area_RNS
#./exe_RNS_CF_peos
rm output_RNS_CF_peos
rm rns_parameter.las
(./exe_RNS_CF_peos > output_RNS_CF_peos ; \
 tail -n 50 rnsphyseq.dat > rns_parameter.las ) &
#(./exe_RNS_WL_peos > output_RNS_WL_peos ; \
# tail -n 50 rnsphyseq.dat > rns_parameter.las ) &
