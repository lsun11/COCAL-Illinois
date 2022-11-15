# set terminal postscript enhanced eps color dashed "Helvetica" 16
# set output 'sotsu3_larg_plot.eps'

set term x11 0
set grid
set title "psi" font "Times-Roman,40"
set xlabel "{/Helvetica=38 x}"
set ylabel "{/Symbol=32 y}"
set mxtics 5
set mytics 2
set xrange [-10:10]
set yrange [-0:3.0]
plot "plot_x.dat" u 1:2 w lp lw 2, \
     "plot_x.dat" u 1:3 w lp lw 2

set term x11 1
set grid
set title "alph" font "Times-Roman,40"
set xlabel "{/Helvetica=38 x}"
set ylabel "{/Symbol=32 y}"
set mxtics 5
set mytics 2
set xrange [-10:10]
set yrange [-0:1.0]
plot "plot_x.dat" u 1:4 w lp lw 2, \
     "plot_x.dat" u 1:5 w lp lw 2

pause -1

