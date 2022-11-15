rm -rf Cocal_debug
cp -rfp ./Cocal Cocal_debug
cd Cocal_debug

cp -rfp work_area_BNS work_area_BBH
cd compile_scripts
sh compile_BBH_CF_circ.sh
#sh compile_BBH_CF_plot_type1.sh

cd ..
cp executable_files/exe_BBH_CF      work_area_BBH/exe_bh_test
#cp executable_files/exe_BBH_CF_plot work_area_BBH/exe_plot
cp ctrl_area_BBH/* work_area_BBH/.

cp -rfp work_area_BBH work_area_BBH_D3

cd work_area_BBH_D3
cp -rfp ../../parameter_sample/bbh_eqm/D3/D* ./

cp exe_bh_test ./D06
cp exe_bh_test ./D07
cp exe_bh_test ./D08
cp exe_bh_test ./D09
cp exe_bh_test ./D10
cp exe_bh_test ./D11
cp exe_bh_test ./D12
cp exe_bh_test ./D13
cp exe_bh_test ./D14
cp exe_bh_test ./D15
cp exe_bh_test ./D16
cp exe_bh_test ./D17
cp exe_bh_test ./D18
cp exe_bh_test ./D19
cp exe_bh_test ./D20
cp exe_bh_test ./D21
cp exe_bh_test ./D22
cp exe_bh_test ./D23
cp exe_bh_test ./D24
cp exe_bh_test ./D25

cd D06
./exe_bh_test >& out_circ_D06 & 
cd ../D07
./exe_bh_test >& out_circ_D07 & 
cd ../D08
./exe_bh_test >& out_circ_D08 & 
cd ../D09
./exe_bh_test >& out_circ_D09 & 
cd ../D10
./exe_bh_test >& out_circ_D10 & 
cd ../D11
./exe_bh_test >& out_circ_D11 &
cd ../D12
./exe_bh_test >& out_circ_D12 & 
cd ../D13
./exe_bh_test >& out_circ_D13 & 
cd ../D14
./exe_bh_test >& out_circ_D14 & 
cd ../D15
./exe_bh_test >& out_circ_D15 & 
cd ../D16
./exe_bh_test >& out_circ_D16 & 
cd ../D17
./exe_bh_test >& out_circ_D17 & 
cd ../D18
./exe_bh_test >& out_circ_D18 & 
cd ../D19
./exe_bh_test >& out_circ_D19 & 
cd ../D20
./exe_bh_test >& out_circ_D20 & 
cd ../D21
./exe_bh_test >& out_circ_D21 & 
cd ../D22
./exe_bh_test >& out_circ_D22 & 
cd ../D23
./exe_bh_test >& out_circ_D23 & 
cd ../D24
./exe_bh_test >& out_circ_D24 & 
cd ../D25
./exe_bh_test >& out_circ_D25 & 

cd ../../..

