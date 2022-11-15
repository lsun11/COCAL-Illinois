set term x11 4

set size sq 1.0,1.0
set zero 1e-20
set pointsize 1.6
set xtics  1
set ytics  1

set xrange [-1.1:1.1]
set yrange [-1.1:1.1] 
#set xrange [-1.0:1.0]
#set yrange [-1.0:1.0] 

plot 'rnsshape.dat' using 1:2 w lp
set autoscale
set size 1, 1

pause -1
