set terminal postscript enhanced eps color dashed "Helvetica" 28 \
fontfile "/usr/share/texmf/fonts/type1/bluesky/cm/cmsy10.pfb"
#set terminal postscript enhanced eps color solid "Helvetica" 16

set output 'mass_radi.eps'
set xlabel "{/Helvetica=24 R [km]}"
#set ylabel "{/Helvetica=28 M}"
set ylabel "{/Helvetica=24 M [M_{/CMSY10 \014}]}"
xx=1
x1 =   9.934201E+00   #  radius of max adm mass
f1 =   2.059835E+00    #  max adm mass

set grid
set size sq
set arrow from x1, graph 0 to x1,f1 nohead lt 2 lw 2 lc rgb "#000000" 
set xrange [8:18]
set yrange [0:2.5]
set mxtics 4 
set mytics 5
set ytics 0,0.5,2.5
set pointsize 1.6
set tics scale  2, 1
set title sprintf("{/Helvetica=24 R=%1.3f     M_{max adm}=%1.3f}", x1,f1)

plot 'ovphy_plot_eosSLy_EOS_00.dat' u 15:11 w l lw 2 notitle,\
     'maxmin_adm_radius.txt' u (($1)*xx):2 w lp  lt 1 lw 4 lc rgb "#000000"  notitle,\
     x<=x1*xx?f1:0/0 lt 2 lw 2 lc rgb "#000000" notitle

