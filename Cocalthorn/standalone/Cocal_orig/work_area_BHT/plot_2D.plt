 set terminal postscript enhanced eps color dashed "Helvetica" 16
 set output 'sotsu2_larg_l=12_2D.eps'
 set xrange [-70:70]
 set yrange [-70:70]
 set zrange [-0.0005:0.0005]
# set view 30,80
 set palette rgbformulae 23,28,3 
 set pm3d map
 splot "rns_contour_xy.dat" using 1:2:3  notitle

