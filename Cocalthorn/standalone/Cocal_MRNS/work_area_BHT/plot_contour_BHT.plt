# Contour plot for conformal factor psi
reset ;
x1=-2;  x2=2;
y1=-2;  y2=2; 

re=0.06248999599680   ;
set grid ;
set size sq;
set xtics x1,1.0,x2 ;
set ytics y1,1.0,y2 ;
set xrange [x1:x2];  
set yrange [y1:y2];
set xlabel "x"  font "Helvetica, 18"
set ylabel "y" font "Helvetica, 18" ; 
set pointsize 0.7 ;
#set size 1.0,1.0 ;
set contour base;
unset surface;
unset key
set view map;
set cntrparam levels incremental  1.0, 0.01, 1.16
#set term postscript portrait enhanced color solid 
set terminal postscript enhanced eps solid color "Helvetica" 24
set output "bht_contour_xy_psi.eps"
#splot "bns_contour_xy_mpt1.dat" u ($1+xns1):2:7 w line palette 
splot "rns_contour_xy_CF.dat" u 1:2:3 w lines


# Contour plot for lapse alph
#reset
#x1=-3;  x2=3;
#y1=-3;  y2=3; xns1=-1.25  ;  xns2=1.25 ;
#re=1.125   ;
#set grid ;
#set size sq;
#set xtics x1,1.0,x2 ;
#set ytics y1,1.0,y2 ;
#set xrange [x1:x2];  set yrange [y1:y2];
#set xlabel "x"  font "Helvetica, 18" ;
#set ylabel "y" font "Helvetica, 18" ;
#set pointsize 0.7 ;
##set size 1.0,1.0 ;
#set contour base;
#unset surface;
#unset key
#set view map;
#set cntrparam levels incremental  0.75, 0.01, 1.0
#set terminal postscript enhanced eps solid color "Helvetica" 24
#set output "bns_contour_xy_mpt1_alph.eps"
#splot "bns_contour_xy_mpt1.dat" u ($1+xns1):2:8 w lines
#
## Contour plot for emd
#reset
#x1=-3;  x2=3;
#y1=-3;  y2=3; xns1=-1.25  ;  xns2=1.25 ;
#re=1.125   ;
#set grid ;
#set size sq;
#set xtics x1,1.0,x2 ;
#set ytics y1,1.0,y2 ;
#set xrange [x1:x2];  set yrange [y1:y2];
#set xlabel "x"   font "Helvetica, 18" ;
#set ylabel "y" font "Helvetica, 18" ;
#set pointsize 0.7 ;
##set size 1.0,1.0 ;
#set contour base;
#unset surface;
#unset key
#set view map;
#set cntrparam levels incremental  0.01, 0.01, 0.09
#set table "ns1.dat"
#splot "bns_contour_xy_mpt1.dat" u ($1+xns1):2:3
#unset table
#set table "ns2.dat"
#splot "bns_contour_xy_mpt2.dat" u ($1+xns2):2:3 
#unset table
#
#reset
#x1=-3;  x2=3.0;
#y1=-3;  y2=3.0 ; xns1=-1.25  ;  xns2=1.25 ;
#rs1=0.625
#rs2=0.625
#
#re=1.125   ;
#set grid ;
#set size sq;
#set xtics x1,1.0,x2 ;
#set ytics y1,1.0,y2 ;
#set xrange [x1:x2];  set yrange [y1:y2];
#set xlabel "x"   font "Helvetica, 18" ;
#set ylabel "y" font "Helvetica, 18" ;
#set pointsize 0.7 ;
##set size 1.0,1.0 ;
#set contour base;
#unset surface;
#unset key
#set view map;
#set cntrparam levels incremental  0.01, 0.01, 0.09
#set terminal postscript enhanced eps solid color "Helvetica" 24
#set output "bns_contour_xy_mpt1_emd.eps"
#plot "bnsshape_xy_mpt1.dat" u (xns1+rs1*($1)):(rs1*($2))  w filledcurve fs solid 0.1 lc rgb '#0000FF' notitle,\
#     "bnsshape_xy_mpt2.dat" u (xns2-rs2*($1)):(rs2*(-($2)))  w filledcurve fs solid 0.1 lc rgb '#0000FF' notitle,\
#     "ns1.dat" w l lc rgb '#FF0000' notitle, "ns2.dat" w l lc rgb '#FF0000' notitle
#
#
## Contour plot for vep
# reset
# x1=-3;  x2=3;
# y1=-3;  y2=3; xns1=-1.25  ;  xns2=1.25 ;
# re=1.125   ;
# set grid ;
# set size sq;
# set xtics x1,1.0,x2 ;
# set ytics y1,1.0,y2 ;
# set xrange [x1:x2];  set yrange [y1:y2];
# set xlabel "x"   font "Helvetica, 18" ;
# set ylabel "y" font "Helvetica, 18" ;
# set pointsize 0.7 ;
# #set size 1.0,1.0 ;
# set contour base;
# unset surface;
# unset key
# set view map;
# set cntrparam levels incremental  -0.1, 0.01, 0.1
# set table "ns3.dat"
# splot "bns_contour_xy_mpt1.dat" u ($1+xns1):2:12
# unset table
# set table "ns4.dat"
# splot "bns_contour_xy_mpt2.dat" u ($1+xns2):2:12
# unset table
#
#reset
#x1=-3;  x2=3.0; 
#y1=-3;  y2=3.0 ; xns1=-1.25  ;  xns2=1.25 ;
#rs1=0.625
#rs2=0.625
#
#re=1.125   ;
#set grid ;
#set size sq;
#set xtics x1,1.0,x2 ;
#set ytics y1,1.0,y2 ;
#set xrange [x1:x2];  set yrange [y1:y2];
#set xlabel "x"   font "Helvetica, 18" ;
#set ylabel "y" font "Helvetica, 18" ;
#set pointsize 0.7 ;
##set size 1.0,1.0 ;
#set contour base;
#unset surface;
#unset key
#set view map;
#set cntrparam levels incremental  -0.1, 0.01, 0.1
#set terminal postscript enhanced eps solid color "Helvetica" 24
#set output "bns_contour_xy_mpt1_vep.eps"
#plot "bnsshape_xy_mpt1.dat" u (xns1+rs1*($1)):(rs1*($2))  w filledcurve fs solid 0.1 lc rgb '#0000FF' notitle,\
#     "bnsshape_xy_mpt2.dat" u (xns2-rs1*($1)):(rs2*(-($2)))  w filledcurve fs solid 0.1 lc rgb '#0000FF' notitle,\
#          "ns3.dat" w l lc rgb '#FF0000' notitle, "ns4.dat" w l lc rgb '#FF0000' notitle
#


#set parametric
#plot "./BBH_contour_xy_new.dat" u 1:2 notitle w lines,  bhp2+re*cos(t),re*sin(t) notitle lt  2, \
#                                           ra1*cos(t),ra1*sin(t) notitle lt -1 lw 1.5, \
#                                           bhp2+ra2*cos(t),ra2*sin(t) notitle lt -1 lw 1.5
