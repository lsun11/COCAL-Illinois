# set terminal postscript enhanced eps color dashed "Helvetica" 16
# set output 'sotsu3_larg_plot.eps'
 set title "Plot of Potential" font "Times-Roman,40"
 set xlabel "{/Helvetica=38 x}"
 set ylabel "{/Symbol=32 y}"
 set grid
 set mxtics 5
 set mytics 2
# set tics font 'Helvetica,22'
 #set yrange [-80:80]
 #set zrange [-0.0005:0.0005]
 #set view 30,80
 #set palette rgbformulae 23,28,3

set logscale x
set logscale y
 set xrange [0.01:1000000]
 set yrange [0.000000001:10]
# set yrange [0.8:2.2]
 set terminal x11 0
 plot "plot_x_mpt1.dat" u 1:(abs($2-1.0)) w lp lw 2, \
      "plot_x_mpt2.dat" u 1:(abs($2-1.0)) w lp lw 2, \
      "plot_x_mpt3.dat" u 1:(abs($2-1.0)) w lp lw 2, 1/x

# set yrange [-0.2:1.2]
 set terminal x11 1
 plot "plot_x_mpt1.dat" u 1:(abs($3-1.0)) w lp lw 2, \
      "plot_x_mpt2.dat" u 1:(abs($3-1.0)) w lp lw 2, \
      "plot_x_mpt3.dat" u 1:(abs($3-1.0)) w lp lw 2, 1/x


# set xrange [-20:20]
# set yrange [0.8:2.2]
 set yrange [0.000000000000001:10]
 set terminal x11 2
 plot "plot_x_mpt1.dat" u 1:(abs($4)) w lp lw 2, \
      "plot_x_mpt2.dat" u 1:(abs($4)) w lp lw 2, \
      "plot_x_mpt3.dat" u 1:(abs($4)) w lp lw 2, 1/x, 1/x**2

# set yrange [-0.2:1.2]
 set terminal x11 3
 plot "plot_x_mpt1.dat" u 1:(abs($5)) w lp lw 2, \
      "plot_x_mpt2.dat" u 1:(abs($5)) w lp lw 2, \
      "plot_x_mpt3.dat" u 1:(abs($5)) w lp lw 2, 1/x, 1/x**2

 set terminal x11 4
 plot "plot_x_mpt1.dat" u 1:(abs($6)) w lp lw 2, \
      "plot_x_mpt2.dat" u 1:(abs($6)) w lp lw 2, \
      "plot_x_mpt3.dat" u 1:(abs($6)) w lp lw 2, 1/x, 1/x**2

pause -1