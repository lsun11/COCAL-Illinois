set terminal postscript enhanced eps color dashed "Helvetica" 24 \
fontfile "/usr/share/texmf/fonts/type1/bluesky/cm/cmsy10.pfb"

set title "iteration100.txt" 
set xlabel "{/Helvetica=28 x}"
set ylabel "{/Symbol=28 y}"
set grid
set mxtics 5
set mytics 5
set xrange [-4:4]
set arrow from  1, graph 0 to  1, graph 1 nohead lt 2 lw 2 lc rgb "#000000"
set arrow from -1, graph 0 to -1, graph 1 nohead lt 2 lw 2 lc rgb "#000000"
set label "star" at -0.3,graph 0.1 rotate by 0 point ps 2
# set yrange [0.8:2.2]
set output 'xaxis_psi_iteration100.eps'
plot "iteration100.txt" u 1:2 w l lw 2 lt 1 lc rgb "red" notitle

set ylabel "{/Symbol=28 a}"
set output 'xaxis_alph_iteration100.eps'
plot "iteration100.txt" u 1:3 w l lw 2 lt 1 lc rgb "red" notitle

set ylabel "{/Symbol=28 b}^{y}"
set output 'xaxis_bvyd_iteration100.eps'
plot "iteration100.txt" u 1:5 w l lw 2 lt 1 lc rgb "red" notitle

set xrange [-2:2]
set ylabel "P/{/Symbol=28 r}"
set output 'xaxis_emd_iteration100.eps'
plot "iteration100.txt" u 1:6 w l lw 2 lt 1 lc rgb "red" notitle
