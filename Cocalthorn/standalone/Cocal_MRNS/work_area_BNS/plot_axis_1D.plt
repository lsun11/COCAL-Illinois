set terminal postscript enhanced eps color dashed "Helvetica" 24 \
fontfile "/usr/share/texmf/fonts/type1/bluesky/cm/cmsy10.pfb"

#set title "iteration5.txt" 
set xlabel "{/Helvetica=28 @^{\\261}r}"
#set ylabel "{/Symbol=28 y}"
set grid
xns1=-1.5
xns2=+1.5
rs1=0.625

set mxtics 5
set mytics 5
set xrange [0:2.5]
set yrange [0.7:1.2]
set arrow from rs1, graph 0 to rs1, graph 1 nohead lt 2 lw 2 lc rgb "#000000"
#set arrow from xns1-1, graph 0 to xns1-1, graph 1 nohead lt 2 lw 2 lc rgb "#000000"
#set arrow from xns2+1, graph 0 to xns2+1, graph 1 nohead lt 2 lw 2 lc rgb "#000000"
#set arrow from xns2-1, graph 0 to xns2-1, graph 1 nohead lt 2 lw 2 lc rgb "#000000"
set label 1 "star" at 1.05*rs1,graph 0.1 rotate by 90 point ps 2
#set label 2 "star" at xns2-0.25,graph 0.1 rotate by 0 point ps 2
# set yrange [0.8:2.2]
set output 'xaxis_alps_1D.eps'
plot "rnsgra_1D_mpt1.ini" u 1:((abs($1)<=rs1)?($2):0/0) w l lw 2 lt 1 lc rgb "#FF0000" title "{/Symbol=28 y}",\
     "rnsgra_1D_mpt1.ini" u 1:2 w l lw 4 lt 0 lc rgb "#FF0000" notitle,\
     "rnsgra_1D_mpt1.ini" u 1:((abs($1)<=rs1)?($3):0/0) w l lw 2 lt 1 lc rgb "#0000FF" title "{/Symbol a}",\
     "rnsgra_1D_mpt1.ini" u 1:3 w l lw 4 lt 0 lc rgb "#0000FF" notitle

set xrange [0:1.2]
set output 'xaxis_emd_1D.eps'
plot "rnsflu_1D_mpt1.ini" u 1:2 w l lw 2 lt 1 lc rgb "#FF0000" title "P/{/Symbol=28 r}_0"
