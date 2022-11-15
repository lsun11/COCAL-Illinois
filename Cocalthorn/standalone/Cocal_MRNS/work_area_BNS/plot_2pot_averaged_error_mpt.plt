# set terminal postscript enhanced eps color dashed "Helvetica" 16
# set output 'sotsu3_larg_plot.eps'
set xlabel "{/Helvetica=38 x}"
set grid

set logscale y
set format y "10^{%L}"
set zero 1e-20
set mxtics 5 ; set mytics 10
set pointsize 1.6
set ticscale  2 1

set xrange [-10:10]
set yrange [0.001:100]
set ylabel "{/Symbol=32 dy}{\Helvetica=32 /}{/Symbol=32 y}{\Helvetica=32 [%]}"
set terminal x11 0
plot "plot_averaged_error_mpt1.dat" u 1:2 w lp lt 0 lw 2, \
     "plot_averaged_error_mpt2.dat" u 1:2 w lp lt 1 lw 2

set xrange [-20:20]
set yrange [0.001:100]
set ylabel "{/Symbol=32 dy}{\Helvetica=32 /}{/Symbol=32 y}{\Helvetica=32 [%]}"
set terminal x11 1
plot "plot_averaged_error_mpt1.dat" u 1:2 w lp lt 0 lw 2, \
     "plot_averaged_error_mpt2.dat" u 1:2 w lp lt 1 lw 2


set xrange [-10:10]
set yrange [0.001:100]
set ylabel "{/Symbol=32 da}{\Helvetica=32 /}{/Symbol=32 a}{\Helvetica=32 [%]}"
set terminal x11 2
plot "plot_averaged_error_mpt1.dat" u 1:3 w lp lt 0 lw 2, \
     "plot_averaged_error_mpt2.dat" u 1:3 w lp lt 1 lw 2

set xrange [-20:20]
set yrange [0.001:100]
set ylabel "{/Symbol=32 da}{\Helvetica=32 /}{/Symbol=32 a}{\Helvetica=32 [%]}"
set terminal x11 3
plot "plot_averaged_error_mpt1.dat" u 1:3 w lp lt 0 lw 2, \
     "plot_averaged_error_mpt2.dat" u 1:3 w lp lt 1 lw 2


pause -1