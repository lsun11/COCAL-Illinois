#set terminal postscript eps color solid 22
#set output 'plot_test1.eps'
set terminal x11 0
set xlabel "\Omega M"
set ylabel "J/M^2"
plot \
     'rnsphyplot_all.dat' u 18:30 w lp

#set output 'plot_test2.eps'
set terminal x11 1
set xlabel "\Omega M"
set ylabel "T/|W|"
plot \
     'rnsphyplot_all.dat' u 18:34 w lp

#set output 'plot_test3.eps'
set terminal x11 2
set xlabel "\Omega M"
set ylabel "M_ADM"
plot \
     'rnsphyplot_all.dat' u 18:20 w lp

pause -1
