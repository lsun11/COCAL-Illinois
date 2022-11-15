set terminal postscript enhanced eps color dashed "Helvetica" 28 \
fontfile "/usr/share/texmf/fonts/type1/bluesky/cm/cmsy10.pfb"
#set terminal postscript enhanced eps color solid "Helvetica" 16

set output 'mass_compa.eps'
set key left bottom
set xlabel "{/Helvetica=28 M/R}"
#set ylabel "{/Helvetica=24 M [M_{/CMSY10 \014}]}"
set ylabel "{/Helvetica=28 M}"
set y2label "{/Helvetica=28 R (Schwarzschild coordinate)}"
xx=1
x1 = 2.146387E-01   # cmpactness at max adm mass
f1 = 2.232491E+00   # max adm mass 
f2 = 1.908993E+00   # adm mass at max compactness
x2 = 2.485391E-01   # max compactness

set grid
set arrow from x1, graph 0 to x1,f1 nohead lt 2 lw 2 lc rgb "#000000"
set arrow from x2, graph 0 to x2,f2 nohead lt 2 lw 2 lc rgb "#000000"
set xrange [0:0.27]
set ytics nomirror 
#set y2range [0.3:1.3]
set mxtics 5 
set mytics 5
set y2tics nomirror
set my2tics 5
set pointsize 1.6
set tics scale  2, 1
set title sprintf("{/Helvetica=28 (M/R)=%1.3f     M_{adm max}=%1.3f}\n {/Helvetica=20 (M/R)_{max}=%1.3f     M_{adm}=%1.3f}", x1,f1,x2,f2)
plot 'ovphy_plot_E2ekadath_EOS_00.dat' u 1:11 w l lt 1 lw 2 lc rgb "#FF0000" axes x1y1 title "{/Helvetica=28 M/R-M}",\
     'ovphy_plot_E2ekadath_EOS_00.dat' u 1:14  w l lt 1 lw 2 lc rgb "#0000FF" axes x1y2 title "{/Helvetica=28 M/R-R}",\
     'maxmin_adm_compa.txt' u 1:2 w lp lt 1 lw 4 lc rgb "#000000" notitle,\
     x<=x1?f1:0/0 lt 2 lw 2 lc rgb "#000000" notitle,\
     'maxmin_compa_adm.txt' u 2:1 w lp lt 1 lw 4 lc rgb "#000000" notitle,\
     x<=x2?f2:0/0 lt 2 lw 2 lc rgb "#000000" notitle


