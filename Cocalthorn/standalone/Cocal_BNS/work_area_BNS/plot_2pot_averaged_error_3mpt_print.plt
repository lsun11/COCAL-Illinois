set terminal postscript enhanced eps color dashed "Helvetica" 24
set xlabel "{/Helvetica=38 x}"
set grid

set logscale y
set format y "10^{%L}"
set zero 1e-20
set mxtics 5 ; set mytics 10
set pointsize 1.6
set ticscale  2 1

set xrange [-15:15]
set yrange [0.001:0.1]
set ylabel "{/Symbol=32 dy}{\Helvetica=32 /}{/Symbol=32 y}{\Helvetica=32 [%]}"
plot "plot_averaged_error_mpt1.dat" u 1:2 w lp lt 0 lw 2, \
     "plot_averaged_error_mpt2.dat" u 1:2 w lp lt 1 lw 2, \
     "plot_averaged_error_mpt3.dat" u 1:2 w lp lt 2 lw 2

set xrange [-15:15]
set yrange [0.001:100]
set ylabel "{/Symbol=32 da}{\Helvetica=32 /}{/Symbol=32 a}{\Helvetica=32 [%]}"
plot "plot_averaged_error_mpt1.dat" u 1:3 w lp lt 0 lw 2, \
     "plot_averaged_error_mpt2.dat" u 1:3 w lp lt 1 lw 2, \
     "plot_averaged_error_mpt3.dat" u 1:3 w lp lt 2 lw 2

pause -1