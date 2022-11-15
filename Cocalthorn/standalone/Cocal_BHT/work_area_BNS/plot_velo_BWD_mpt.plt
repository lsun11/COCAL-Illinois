# Gnuplot script file for plotting fluid velocity in the XY plane
# For the NS shape plot in COCP2 (Xcm,Ycm)=(dis-x2, -y2)
# For the velocity inside COCP2 the vector direction is taken care inside code.
#
reset ;
x1=-3.0;   x2=3.0; 
y1=-3.0;   y2=3.0 ; 
xns1=-1.5  ;  xns2=1.5 ;
re=1.125 ; 
set size sq 
set xtics x1,1.0,x2 ;
set ytics y1,1.0,y2 ;
set grid;  set xrange [x1:x2];  set yrange [y1:y2];
set xlabel "x"  font "Helvetica, 18" ;
set ylabel "y" font "Helvetica, 18" ; 
set pointsize 0.7 ;
set parametric;
mf=40 ;
set terminal postscript enhanced eps solid color "Helvetica" 24
set output "bns_contour_xy_mpt1_vxy.eps"
#plot "./BBH_contour_xy.dat" u 1:2:(mf*($5)):(mf*($6)) w vec notitle lt 1 lw 2, bhp2+re*cos(t),re*sin(t) notitle lt  2, \
#                                                                          ra1*cos(t),ra1*sin(t) notitle lt -1 lw 1.5, \
#                                                                     bhp2+ra2*cos(t),ra2*sin(t) notitle lt -1 lw 1.5
plot "bnsshape_xy_mpt1.dat" u (xns1+$1):2  w filledcurve fs solid 0.1 lc rgb '#0000FF' notitle,\
     "bnsshape_xy_mpt2.dat" u (xns2-$1):(-($2))  w filledcurve fs solid 0.1 lc rgb '#0000FF' notitle,\
     "bns_contour_xy_mpt1.dat" u (xns1+$1):2:(mf*($4)):(mf*($5)) every 2:2 w vec notitle lt 1 lw 2,\
     "bns_contour_xy_mpt2.dat" u (xns2-$1):(-($2)):(mf*($4)):(mf*($5)) every 2:2 w vec notitle lt 1 lw 2
