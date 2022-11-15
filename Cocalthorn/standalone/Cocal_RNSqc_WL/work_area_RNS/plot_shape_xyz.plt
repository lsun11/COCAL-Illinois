set size sq 1.0,1.0
set zero 1e-20
set pointsize 1.6
set xtics  1
set ytics  1

set xrange [-1.1:1.1]
set yrange [-1.1:1.1] 
#set xrange [-1.0:1.0]
#set yrange [-1.0:1.0] 

set term x11 1
plot 'rnsshape_seq_xy.dat' using 1:2 w lp
set term x11 2
plot 'rnsshape_seq_xz.dat' using 1:2 w lp
set term x11 3
plot 'rnsshape_seq_zy.dat' using 1:2 w lp

set autoscale
set size 1, 1

pause -1
