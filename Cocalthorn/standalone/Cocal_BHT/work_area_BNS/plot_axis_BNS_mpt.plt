set terminal postscript enhanced eps color dashed "Helvetica" 24 \
fontfile "/usr/share/texmf/fonts/type1/bluesky/cm/cmsy10.pfb"

set title "iteration85.txt" 
set xlabel "{/Helvetica=28 x}"
set ylabel "{/Symbol=28 y}"
set grid
xns1=-1.25
xns2=+1.25
rs1=0.625
rs2=0.625

set mxtics 5
set mytics 5
set xrange [-4:4]
set arrow from xns1+rs1, graph 0 to xns1+rs1, graph 1 nohead lt 2 lw 2 lc rgb "#000000"
set arrow from xns1-rs1, graph 0 to xns1-rs1, graph 1 nohead lt 2 lw 2 lc rgb "#000000"
set arrow from xns2+rs2, graph 0 to xns2+rs2, graph 1 nohead lt 2 lw 2 lc rgb "#000000"
set arrow from xns2-rs2, graph 0 to xns2-rs2, graph 1 nohead lt 2 lw 2 lc rgb "#000000"
set label 1 "star" at xns1-0.25,graph 0.1 rotate by 0 point ps 2
set label 2 "star" at xns2-0.25,graph 0.1 rotate by 0 point ps 2
# set yrange [0.8:2.2]
set output 'xaxis_psi_iteration85.eps'
plot "iteration85_mpt1.txt" u 1:2 w l lw 2 lt 1 lc rgb "#FF0000" notitle,\
     "iteration85_mpt2.txt" u 1:2 w l lw 2 lt 1 lc rgb "#0000FF" notitle,\
     "iteration85_mpt3.txt" u 1:2 w l lw 2 lt 1 lc rgb "#008000" notitle

set ylabel "{/Symbol=28 a}"
unset label 1
unset label 2
set label 1 "star" at xns1-0.25,graph 0.9 rotate by 0 point ps 2
set label 2 "star" at xns2-0.25,graph 0.9 rotate by 0 point ps 2
set output 'xaxis_alph_iteration85.eps'
plot "iteration85_mpt1.txt" u 1:3 w l lw 2 lt 1 lc rgb "#FF0000" notitle,\
     "iteration85_mpt2.txt" u 1:3 w l lw 2 lt 1 lc rgb "#0000FF" notitle,\
     "iteration85_mpt3.txt" u 1:3 w l lw 2 lt 1 lc rgb "#008000" notitle


set ylabel "{/Symbol=28 b}^{y}"
set output 'xaxis_bvyd_iteration85.eps'
plot "iteration85_mpt1.txt" u 1:5 w l lw 2 lt 1 lc rgb "#FF0000" notitle,\
     "iteration85_mpt2.txt" u 1:5 w l lw 2 lt 1 lc rgb "#0000FF" notitle,\
     "iteration85_mpt3.txt" u 1:5 w l lw 2 lt 1 lc rgb "#008000" notitle


set ylabel "P/{/Symbol=28 r}"
unset label 1
unset label 2
set label 1 "star" at xns1-0.25,graph 0.1 rotate by 0 point ps 2
set label 2 "star" at xns2-0.25,graph 0.1 rotate by 0 point ps 2
set output 'xaxis_emd_iteration85.eps'
plot "iteration85_mpt1.txt" u 1:((abs($1-xns1)<=rs1)?($6):0/0) w l lw 2 lt 1 lc rgb "#FF0000" notitle,\
     "iteration85_mpt2.txt" u 1:((abs($1-xns2)<=rs2)?($6):0/0) w l lw 2 lt 1 lc rgb "#0000FF" notitle

set ylabel "vep"
unset label 1
unset label 2
set label 1 "star" at xns1-0.25,graph 0.1 rotate by 0 point ps 2
set label 2 "star" at xns2-0.25,graph 0.1 rotate by 0 point ps 2
set output 'xaxis_vep_iteration85.eps'
plot "iteration85_mpt1.txt" u 1:((abs($1-xns1)<=rs1)?($7):0/0) w l lw 2 lt 1 lc rgb "#FF0000" notitle,\
     "iteration85_mpt2.txt" u 1:((abs($1-xns2)<=rs2)?($7):0/0) w l lw 2 lt 1 lc rgb "#0000FF" notitle
