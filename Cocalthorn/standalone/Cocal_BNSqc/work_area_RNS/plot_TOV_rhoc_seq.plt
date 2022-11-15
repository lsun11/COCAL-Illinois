set terminal postscript enhanced eps color dashed "Helvetica" 24 \
fontfile "/usr/share/texmf/fonts/type1/bluesky/cm/cmsy10.pfb"
#set terminal postscript enhanced eps color solid "Helvetica" 16

set output 'mass_rhoc.eps'
set key right bottom
set xlabel "{/Helvetica=24 rest mass density {/Symbol r}_c }"
#set xlabel "{/Helvetica=24 rest mass density log{/Symbol r}_c }"
#set ylabel "{/Helvetica=28 M }"
set ylabel "{/Helvetica=24 M [M_{/CMSY10 \014}]}"
#
# Run findstats.f90 as follows
# ./a.out ovphy_plot_K123.6_EOS_00.dat 3 11 => cut last line in maxmin_adm_rho0.txt
# modify x2,f2
# ./a.out ovphy_plot_K123.6_EOS_00.dat 3  9 => cut last line in maxmin_restmass_rho0.txt
# modify x1,f1
#
x1 = 1.590000E+15   # rest mass density xmax
f1 = 1.820253E+00   # admmass at xmax

# Scale x-axis (density) by scale factor sc.
sc   = 1.0E+15
x1sc = x1/sc
xa   = 0.01*x1/sc
xb   = 1.5*x1/sc
set xlabel "{/Helvetica=24 rest mass density {/Symbol r}_c ({/Symbol \264} 10^{15})}"

set grid
set size sq 
#set logscale x
#set format x "10^{%L}"
#set xrange [x1/100:10*x1]
set xrange [xa:xb]
set yrange [0:1.3*f1]
#set ytics 0,0.5,2.5
set arrow 1 from x1sc, graph 0 to x1sc,f1 nohead lt 2 lw 1 lc rgb "#000000" 
#set ytics 0.5
set mxtics 5 
set mytics 5
set pointsize 1.6
set tics scale  2, 1
#set title sprintf("{/Helvetica=24 {/Symbol r}_{max}=%9.3e     M_{max adm}=%1.3f}", x1,f1)
        
plot 'ovphy_plot_K123.6_EOS_00.dat' u (($4)/sc):11 w l lt 1 lw 2 lc rgb "#FF0000"  notitle,\
     'maxmin_adm_rho0.txt' u (($1)/sc):2 w lp  lt 1 lw 4 lc rgb "#000000"  notitle,\
     x<=(x1sc)?f1:0/0 lt 2 lw 1 lc rgb "#000000" notitle


