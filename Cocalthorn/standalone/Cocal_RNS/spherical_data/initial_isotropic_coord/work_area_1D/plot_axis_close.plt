set format x "10^{%L}"
set zero 1e-20
set mxtics 10 ; set mytics 4
set pointsize 1.6
set ticscale  2 1

call "plot_XYZ_close.plt" "frggraplot.dat" "2" "x/R" "psi" 
call "plot_XYZ_close.plt" "frggraplot.dat" "3" "x/R" "alpha" 
call "plot_XYZ_close.plt" "frggraplot.dat" "4" "x/R" "p/rho" 

pause -1
