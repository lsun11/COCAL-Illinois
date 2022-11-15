 set terminal postscript enhanced eps color dashed "Helvetica" 16
 set output 'fig_wave_form.eps'
 set xrange [-80:80]
 set yrange [-80:80]
 set zrange [-0.0001:0.0001]
 set view 30,80
 set palette rgbformulae 23,28,3
 a = 0.167
 splot "rns_contour_xy_mpt3.dat" using 1:2:(($3)+a/sqrt(($1)**2+($2)**2)) with pm3d notitle

