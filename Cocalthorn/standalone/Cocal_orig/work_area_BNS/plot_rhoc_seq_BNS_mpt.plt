set terminal postscript enhanced eps color dashed "Helvetica" 24 \
fontfile "/usr/share/texmf/fonts/type1/bluesky/cm/cmsy10.pfb"

set ylabel "{/Helvetica=28 M_0}"
set xlabel "{/Helvetica=28 {/Symbol=28 r}_c}"
set grid
x1 = 3.240065e-01   # restmass xmax
f1 = 1.797995e-01   # restmass at xmax
x2 = 3.240065E-01   # corot restmass xmax
f2 = 1.869095E-01   # corot restmass at xmax
set mxtics 5
set mytics 5
set arrow from x1, graph 0 to x1,f1 nohead lt 2 lw 2 lc rgb "#000000"
set arrow from x2, graph 0 to x2,f2 nohead lt 2 lw 2 lc rgb "#000000"
set xrange [0:0.5]
set yrange [0:0.24]
set output 'restmass_rhoc_BNS.eps'
plot "ovphy_plot_TOV_seq_EOS_00.dat" u 3:9  w l  lt 1 lw 2 lc rgb "#0000FF" title "TOV",\
     "maxmin_mrest_rhoc.txt"         u 1:2  w lp lt 1 lw 4 lc rgb "#000000" notitle,\
     "sequence_data.txt"             u 11:3 w l lw 2 lt 1 lc rgb "#FF0000" notitle,\
     "maxmin_mrest_tov_rho0.txt"     u 1:2  w lp lt 1 lw 4 lc rgb "#000000" notitle,\
     x<=x1?f1:0/0 lt 2 lw 2 lc rgb "#000000" notitle,\
     x<=x2?f2:0/0 lt 2 lw 2 lc rgb "#000000" notitle


reset
#set ylabel "{/Helvetica=28 J}"
#set xlabel "{/Helvetica=28 M_0{/Symbol=28 W}}"
set ylabel "{/Helvetica=28 M_0}"
set xlabel "{/Helvetica=28 M/R}"
set grid
set mxtics 5
set mytics 5
#set xrange [0:0.6]
#set yrange [0:0.22]
set output 'compa_mrest_BNS.eps'
plot "sequence_data.txt" u 10:3 w l lw 2 lt 1 lc rgb "#FF0000" notitle

