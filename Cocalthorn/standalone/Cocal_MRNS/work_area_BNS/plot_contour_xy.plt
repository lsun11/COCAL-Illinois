# Gnuplot script file for plotting contours of Pse and LPse in the XY plane

reset ;
x1=-3.0;  x2=5.0; 
y1=-3.0;  y2=3.0 ; 

bhp1=0.0  ;  bhp2=2.5 ;
re=1.125   ;  
ra1=0.25   ;  ra2=0.25   ;
#rb1=1.25  ;  rb2=1.25  ;

set xtics x1,1.0,x2 ;
set ytics y1,1.0,y2 ;
#set mytics 10 ;
#set format y "10^{%L}" ; 
set grid;  set xrange [x1:x2];  set yrange [y1:y2];
#set size 0.8,0.5 ;

set xlabel "x"  0,0 font "Helvetica, 18" ;
set ylabel "y" font "Helvetica, 18" ; 

set pointsize 0.7 ;

set contour base;
unset surface;
#set contour
set view 0,0;

#Data for the LPse contours
#set cntrparam levels  incremental 0.4, 0.05, 1.0

#Data for the Pse contours
set cntrparam levels incremental  1.0, 0.1, 2.9

#set term table
set term postscript portrait enhanced color solid 
set output "psi_contour.ps"

#set output "contour.txt"
#For plotting Pse we need 1:2:3 and for LPse 1:2:4 
splot "./BBH_contour_xy_new.dat" u 1:2:3 w lines  

#set term x11
#set term postscript portrait enhanced color solid 
#set output "psi_contour.ps"

#set parametric
#plot "./BBH_contour_xy_new.dat" u 1:2 notitle w lines,  bhp2+re*cos(t),re*sin(t) notitle lt  2, \
#                                           ra1*cos(t),ra1*sin(t) notitle lt -1 lw 1.5, \
#                                           bhp2+ra2*cos(t),ra2*sin(t) notitle lt -1 lw 1.5
