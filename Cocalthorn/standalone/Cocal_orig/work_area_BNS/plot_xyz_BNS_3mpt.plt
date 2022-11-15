set terminal postscript enhanced eps color dashed "Helvetica" 24 \
fontfile "/usr/share/texmf/fonts/type1/bluesky/cm/cmsy10.pfb"

# set title "Plot of Potential" font "Times-Roman,40"
 xns1=-1.25
 xns2=+1.25
 rs1=0.6
 rs2=0.6

 set xlabel "{/Helvetica=38 x}"
 set grid
 set mxtics 5
 set mytics 2
# set tics font 'Helvetica,22'
 #set yrange [-80:80]
 #set zrange [-0.0005:0.0005]
 #set view 30,80
 #set palette rgbformulae 23,28,3

 set xrange [-10:10]
# set yrange [0.8:2.2]
 set output 'psi_x_plot.eps'
 set ylabel "{/Symbol=32 y}"
 plot "plot_x_mpt1.dat" u 1:2 w l lw 2 notitle, \
      "plot_x_mpt2.dat" u 1:2 w l lw 2 notitle, \
      "plot_x_mpt3.dat" u 1:2 w l lw 2 notitle

# set yrange [-0.2:1.2]
 set output 'alph_x_plot.eps'
 set ylabel "{/Symbol=32 a}"
 plot "plot_x_mpt1.dat" u 1:3 w l lw 2 notitle, \
      "plot_x_mpt2.dat" u 1:3 w l lw 2 notitle, \
      "plot_x_mpt3.dat" u 1:3 w l lw 2 notitle


# set xrange [-20:20]
# set yrange [0.8:2.2]
 set output 'bvxd_x_plot.eps'
 set ylabel "{/Symbol=32 b}^x"
 plot "plot_x_mpt1.dat" u 1:4 w l lw 2 notitle, \
      "plot_x_mpt2.dat" u 1:4 w l lw 2 notitle, \
      "plot_x_mpt3.dat" u 1:4 w l lw 2 notitle

# set yrange [-0.2:1.2]
 set output 'bvyd_x_plot.eps'
 set ylabel "{/Symbol=32 b}^y"
 plot "plot_x_mpt1.dat" u 1:5 w l lw 2 notitle, \
      "plot_x_mpt2.dat" u 1:5 w l lw 2 notitle, \
      "plot_x_mpt3.dat" u 1:5 w l lw 2 notitle

 set output 'bvzd_x_plot.eps'
 set ylabel "{/Symbol=32 b}^z"
 plot "plot_x_mpt1.dat" u 1:6 w l lw 2 notitle, \
      "plot_x_mpt2.dat" u 1:6 w l lw 2 notitle, \
      "plot_x_mpt3.dat" u 1:6 w l lw 2 notitle

 set xrange [-2.5:2.5]
 set yrange [0:0.12]
 set output 'emd_x_plot.eps'
 set ylabel "{/Symbol=32 r}"
 plot "plot_x_mpt1.dat" u 1:((abs($1-xns1)<=rs1)?($7):0/0) w l lw 2 notitle, \
      "plot_x_mpt2.dat" u 1:((abs($1-xns2)<=rs2)?($7):0/0) w l lw 2 notitle
