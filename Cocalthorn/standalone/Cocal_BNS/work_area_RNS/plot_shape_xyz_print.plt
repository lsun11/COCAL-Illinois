set terminal postscript enhanced eps color dashed "Helvetica" 24
#set term postscript eps color solid 22
set size sq 1.0,1.0
set zero 1e-20
set pointsize 1.6
set xtics  1
set ytics  1
set grid
set xrange [-1.2:1.2]
set yrange [-1.2:1.2]
#set xrange [-1.0:1.0]
#set yrange [-1.0:1.0] 

set output "rnsshape_seq_xy.eps"
plot 'rnsshape_seq_xy.dat' using 1:2 every :::0::6 w lp lw 1 lt 1 lc rgb "#FF0000" notitle

set output "rnsshape_seq_xz.eps"
plot 'rnsshape_seq_xz.dat' using 1:2 every :::0::6 w lp lw 1 lt 1 lc rgb "#FF0000" notitle

set output "rnsshape_seq_zy.eps"
plot 'rnsshape_seq_zy.dat' using 2:1 every :::0::6 w lp lw 1 lt 1 lc rgb "#FF0000" notitle


#every A:B:C:D:E:F
#
#A: line increment
#B: data block increment
#C: The first line
#D: The first data block
#E: The last line
#F: The last data block
