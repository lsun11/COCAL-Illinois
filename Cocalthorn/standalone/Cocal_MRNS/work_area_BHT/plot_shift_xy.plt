# Gnuplot script file for plotting vector plot in the XY plane

reset ;
x1=-3;  x2=5.0; 
y1=-3;   y2=3.0 ; 

bhp1=0.0  ;  bhp2=2.5 ;
re=1.125 ;  
ra1=0.25   ;  ra2=0.25   ;
#rb1=1.25  ;  rb2=1.25  ;

set xtics x1,0.5,x2 ;
set ytics y1,0.5,y2 ;
#set mytics 10 ;
set grid;  set xrange [x1-0.5:x2+0.5];  set yrange [y1-0.5:y2+0.5];
#set size 0.8,0.5 ;

    
set xlabel "x"  0,0 font "Helvetica, 18" ;
set ylabel "y" font "Helvetica, 18" ; 

set pointsize 0.7 ;

set parametric;
mf=5 ;

#set term postscript portrait enhanced color solid 
#set output "bvxy.ps"

plot "./BBH_contour_xy.dat" u 1:2:(mf*($5)):(mf*($6)) w vec notitle lt 1 lw 2, bhp2+re*cos(t),re*sin(t) notitle lt  2, \
                                                                          ra1*cos(t),ra1*sin(t) notitle lt -1 lw 1.5, \
                                                                     bhp2+ra2*cos(t),ra2*sin(t) notitle lt -1 lw 1.5
