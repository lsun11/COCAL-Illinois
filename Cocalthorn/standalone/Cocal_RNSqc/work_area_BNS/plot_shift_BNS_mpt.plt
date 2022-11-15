# Gnuplot script file for plotting shift in the XY plane wrt the Central CS (ARCP)
# For the NS shape plot in COCP2 (Xcm,Ycm)=(dis-x2, -y2)
# For the shift inside COCP2 the vector direction is taken care inside code.

reset ;
x1=-3;   x2=3.0; 
y1=-3;   y2=3.0 ; 
xns1=-1.25  ;  xns2=1.25 ;
rs1=0.625
rs2=0.625

re=1.125 ;  
set size sq;
set xtics x1,1.0,x2 ;
set ytics y1,1.0,y2 ;
set grid;  set xrange [x1:x2];  set yrange [y1:y2];
set xlabel "x"  font "Helvetica, 18" ;
set ylabel "y" font "Helvetica, 18" ; 
set pointsize 0.7 ;
set parametric;
mf=8 ;
set style arrow 1 head filled size screen 0.01,15,30 lt 1 lw 2
set terminal postscript enhanced eps solid color "Helvetica" 24
set output "bns_contour_xy_mpt1_bvxy.eps"
#plot "./BBH_contour_xy.dat" u 1:2:(mf*($5)):(mf*($6)) w vec notitle lt 1 lw 2, bhp2+re*cos(t),re*sin(t) notitle lt  2, \
#                                                                          ra1*cos(t),ra1*sin(t) notitle lt -1 lw 1.5, \
#                                                                     bhp2+ra2*cos(t),ra2*sin(t) notitle lt -1 lw 1.5
plot "bnsshape_xy_mpt1.dat" u (xns1+rs1*($1)):(rs1*($2))  w filledcurve fs solid 0.1 lc rgb '#0000FF' notitle,\
     "bnsshape_xy_mpt2.dat" u (xns2-rs1*($1)):(rs2*(-($2)))  w filledcurve fs solid 0.1 lc rgb '#0000FF' notitle,\
     "bns_contour_xy_mpt1.dat" u ($1+xns1):2:(mf*($9)):(mf*($10)) every 6:6 w vec notitle arrowstyle 1
