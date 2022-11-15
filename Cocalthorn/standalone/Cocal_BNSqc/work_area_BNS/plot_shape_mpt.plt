set terminal postscript enhanced eps color dashed "Helvetica" 24 \
fontfile "/usr/share/texmf/fonts/type1/bluesky/cm/cmsy10.pfb"

set grid
set size sq 1.0,1.0
set zero 1e-20
set pointsize 1.6
set xtics  1
set ytics  1

set xrange [-1.2:1.2]
set yrange [-1.2:1.2]
#set xrange [-1.0:1.0]
#set yrange [-1.0:1.0]

set output 'xy_shape.eps'
plot 'rnsshape_mpt1.dat' using 1:2 every :::0::0 w lp title "xy shape"

set output 'xz_shape.eps'
plot 'rnsshape_mpt1.dat' using 1:2 every :::1::1 w lp title "xz shape"

set output 'yz_shape.eps'
plot 'rnsshape_mpt1.dat' using 1:2 every :::2::2 w lp title "yz shape"

