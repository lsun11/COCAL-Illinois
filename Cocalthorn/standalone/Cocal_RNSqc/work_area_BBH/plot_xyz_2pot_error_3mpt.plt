# set terminal postscript enhanced eps color dashed "Helvetica" 16
# set output 'sotsu3_larg_plot.eps'
set xlabel "{/Helvetica=38 x}"
set ylabel "{/Symbol=32 y}"
set grid
set mxtics 5
set mytics 2

set logscale y
set format y "10^{%L}"
set zero 1e-20
set mxtics 4 ; set mytics 10
set pointsize 1.6
set ticscale  2 1

set xrange [-10:10]
set yrange [0.001:100]
set terminal x11 0
plot "plot_x_mpt1.dat" u 1:6 w lp lt 0 lw 2, \
     "plot_x_mpt2.dat" u 1:6 w lp lt 1 lw 2, \
     "plot_x_mpt3.dat" u 1:6 w lp lt 2 lw 2

set xrange [-20:20]
set yrange [0.001:100]
set terminal x11 1
plot "plot_x_mpt1.dat" u 1:6 w lp lt 0 lw 2, \
     "plot_x_mpt2.dat" u 1:6 w lp lt 1 lw 2, \
     "plot_x_mpt3.dat" u 1:6 w lp lt 2 lw 2


set xrange [-10:10]
set yrange [0.001:100]
set terminal x11 2
plot "plot_x_mpt1.dat" u 1:7 w lp lt 0 lw 2, \
     "plot_x_mpt2.dat" u 1:7 w lp lt 1 lw 2, \
     "plot_x_mpt3.dat" u 1:7 w lp lt 2 lw 2

set xrange [-20:20]
set yrange [0.001:100]
set terminal x11 3
plot "plot_x_mpt1.dat" u 1:7 w lp lt 0 lw 2, \
     "plot_x_mpt2.dat" u 1:7 w lp lt 1 lw 2, \
     "plot_x_mpt3.dat" u 1:7 w lp lt 2 lw 2


pause -1