# set terminal postscript enhanced eps color dashed "Helvetica" 16
# set output 'sotsu3_larg_plot.eps'
 set title "Plot of Potential" font "Times-Roman,40"
 set xlabel "{/Helvetica=38 x}"
 set ylabel "{/Symbol=32 y}"
 set mxtics 5
 set mytics 2
# set tics font 'Helvetica,22'
 set xrange [-10:10]
 set yrange [-0:2.2]
 #set yrange [-80:80]
 #set zrange [-0.0005:0.0005]
 #set view 30,80
 #set palette rgbformulae 23,28,3
 plot "plot_x.dat" u 1:2 w lp lw 2, \
      "plot_x.dat" u 1:3 w lp lw 2

pause -1