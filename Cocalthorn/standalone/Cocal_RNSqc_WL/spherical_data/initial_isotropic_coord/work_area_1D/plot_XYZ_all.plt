set term x11 $1
#set terminal postscript eps color solid 22
#set output 'metric_flus.eps'
#set linestyle 20 lt -1 lw 1
#set key box linestyle 20
set logscale x
set xlabel "$2"
set ylabel "$3"
plot '$0' u 1:$1  title '$3' w l
set nologscale x
#pause -1
