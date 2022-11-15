LD_FLAGS=-L/home/antonios/repositories/cocal_illinois/lib
#LD_FLAGS=-L/u/sciteam/tsokaros/COCAL_ILLINOIS/lib
#LD_FLAGS=-L/u/sciteam/ruiz1/COCAL_ILLINOIS/Antonios_lastversion/cocal_illinois/lib

#all: read_bin_mrns
#
#read_bin_nsns: read_Cocal_nsns.C
#	c++ -o read_bin_nsns read_Cocal_nsns.C $(LD_FLAGS) -lcocalnsns -lgfortran
#
#read_bin_rns: read_Cocal_rns.C
#	c++ -o read_bin_rns read_Cocal_rns.C $(LD_FLAGS) -lcocalrns -lgfortran

all: read_bin_nsns  read_bin_nsns_pw  read_bin_nsns_pw_qc  read_bin_rns  read_bin_rns_pw  read_bin_rns_pw_qc  read_bin_rns_pw_qc_wl  read_bin_bht_pw  read_bin_mrns


read_bin_nsns : read_Cocal_nsns.C
	c++ -o read_bin_nsns read_Cocal_nsns.C $(LD_FLAGS) -lcocalnsns -lgfortran

read_bin_nsns_pw : read_Cocal_nsns_pw.c
	c++ -o read_bin_nsns_pw read_Cocal_nsns_pw.c $(LD_FLAGS) -lcocalnsns -lgfortran

read_bin_nsns_pw_qc : read_Cocal_nsns_pw_qc.c
	c++ -o read_bin_nsns_pw_qc read_Cocal_nsns_pw_qc.c $(LD_FLAGS) -lcocalnsnsqc -lgfortran


read_bin_rns : read_Cocal_rns.C
	c++ -o read_bin_rns read_Cocal_rns.C $(LD_FLAGS) -lcocalrns -lgfortran

read_bin_rns_pw : read_Cocal_rns_pw.C
	c++ -o read_bin_rns_pw  read_Cocal_rns_pw.C $(LD_FLAGS) -lcocalrns -lgfortran

read_bin_rns_pw_qc : read_Cocal_rns_pw_qc.C
	c++ -o read_bin_rns_pw_qc  read_Cocal_rns_pw_qc.C $(LD_FLAGS) -lcocalrnsqc -lgfortran

read_bin_rns_pw_qc_wl : read_Cocal_rns_pw_qc_WL.C
	c++ -o read_bin_rns_pw_qc_wl  read_Cocal_rns_pw_qc_WL.C $(LD_FLAGS) -lcocalrnsqc_wl -lgfortran


read_bin_bht_pw : read_Cocal_bht_pw.C
	c++ -o read_bin_bht_pw  read_Cocal_bht_pw.C $(LD_FLAGS) -lcocalbht -lgfortran


read_bin_mrns : read_Cocal_mrns.C
	c++ -o read_bin_mrns read_Cocal_mrns.C $(LD_FLAGS) -lcocalmrns -lgfortran

