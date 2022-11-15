set terminal postscript enhanced eps color dashed "Helvetica" 24 \
fontfile "/usr/share/texmf/fonts/type1/bluesky/cm/cmsy10.pfb"

#set title "iteration5.txt" 
#set xlabel "{/Helvetica=28 @^{\\261}r}"
set xlabel "{/Helvetica=28 x}"
set ylabel "{/Helvetica=28 HaC}"
set grid
rs1=0.7597667
xcm=1.25

#set mxtics 5
#set mytics 5
set xrange [-xcm:xcm]
#set yrange [0.7:1.2]
set arrow from rs1, graph 0 to rs1, graph 1 nohead lt 2 lw 2 lc rgb "#000000"
set arrow from -rs1, graph 0 to -rs1, graph 1 nohead lt 2 lw 2 lc rgb "#000000"
set arrow from xcm, graph 0 to xcm, graph 1 nohead lt 2 lw 2 lc rgb "#000000"
#set label 1 "star" at 1.05*rs1,graph 0.1 rotate by 90 point ps 2
#set label 2 "star" at xns2-0.25,graph 0.1 rotate by 0 point ps 2
set output 'Hac_vio_1D.eps'
plot "HaC_phi000_mpt1.txt" u 1:(log10(abs($2))) w lp lw 1 lt 1 lc rgb "#FF0000" title "HaC",\
     "HaC_phi180_mpt1.txt" u (-($1)):(log10(abs($2))) w lp lw 1 lt 1 lc rgb "#FF0000" notitle


set ylabel "{/Helvetica=28 MoC}"
set output 'MoC_vio_1D.eps'
plot "MoC_by_phi000_mpt1.txt" u 1:(log10(abs($2))) w lp lw 1 lt 1 lc rgb "#FF0000" title "{/Symbol=28 b}^y",\
     "MoC_by_phi180_mpt1.txt" u (-($1)):(log10(abs($2))) w lp lw 1 lt 1 lc rgb "#FF0000" notitle 
