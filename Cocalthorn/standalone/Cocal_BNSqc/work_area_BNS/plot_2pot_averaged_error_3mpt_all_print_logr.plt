set terminal postscript enhanced eps color dashed "Helvetica" 24
#set terminal postscript enhanced eps color solid "Helvetica" 16

set xlabel "{/Helvetica=38 r}"
set grid

set logscale x
set logscale y
set format x "10^{%L}"
set format y "10^{%L}"
set zero 1e-20
set mxtics 10 ; set mytics 10
set pointsize 1.6
set tics scale  2, 1

set xrange [0.1:1000000]

set yrange [0.0001:10]
set ylabel "{/Helvetica=44 <}{/Helvetica=32 |}\
{/Symbol=32 dy}{/Helvetica=32 /}{/Symbol=32 y}\
{/Helvetica=32 |}{/Helvetica=44 > }{/Helvetica=32 [%]}"
set output 'plot_2pot_averaged_error_3mpt_psi_logr.eps'
plot \
"../work_area_poisson_test_G1/plot_averaged_error_mpt1.dat" \
u ( 1.25+($1)):2 w l lt 1 lw 2 notitle, \
"../work_area_poisson_test_G1/plot_averaged_error_mpt2.dat" \
u (-1.25+($1)):2 w l lt 1 lw 2 notitle, \
"../work_area_poisson_test_G1/plot_averaged_error_mpt3.dat" \
u ( 0.00+($1)):2 w l lt 1 lw 2 notitle, \
"../work_area_poisson_test_G2/plot_averaged_error_mpt1.dat" \
u ( 1.25+($1)):2 w l lt 2 lw 2 notitle, \
"../work_area_poisson_test_G2/plot_averaged_error_mpt2.dat" \
u (-1.25+($1)):2 w l lt 2 lw 2 notitle, \
"../work_area_poisson_test_G2/plot_averaged_error_mpt3.dat" \
u ( 0.00+($1)):2 w l lt 2 lw 2 notitle, \
"../work_area_poisson_test_G3/plot_averaged_error_mpt1.dat" \
u ( 1.25+($1)):2 w l lt 3 lw 2 notitle, \
"../work_area_poisson_test_G3/plot_averaged_error_mpt2.dat" \
u (-1.25+($1)):2 w l lt 3 lw 2 notitle, \
"../work_area_poisson_test_G3/plot_averaged_error_mpt3.dat" \
u ( 0.00+($1)):2 w l lt 3 lw 2 notitle, \
"../work_area_poisson_test_G4/plot_averaged_error_mpt1.dat" \
u ( 1.25+($1)):2 w l lt 4 lw 2 notitle, \
"../work_area_poisson_test_G4/plot_averaged_error_mpt2.dat" \
u (-1.25+($1)):2 w l lt 4 lw 2 notitle, \
"../work_area_poisson_test_G4/plot_averaged_error_mpt3.dat" \
u ( 0.00+($1)):2 w l lt 4 lw 2 notitle


set yrange [0.0001:100]
set ylabel "{/Helvetica=44 <}{/Helvetica=32 |}\
{/Symbol=32 da}{/Helvetica=32 /}{/Symbol=32 a}\
{/Helvetica=32 |}{/Helvetica=44 > }{/Helvetica=32 [%]}"
set output 'plot_2pot_averaged_error_3mpt_alpha_logr.eps'
plot \
"../work_area_poisson_test_G1/plot_averaged_error_mpt1.dat" \
u ( 1.25+($1)):3 w l lt 1 lw 2 notitle, \
"../work_area_poisson_test_G1/plot_averaged_error_mpt2.dat" \
u (-1.25+($1)):3 w l lt 1 lw 2 notitle, \
"../work_area_poisson_test_G1/plot_averaged_error_mpt3.dat" \
u ( 0.00+($1)):3 w l lt 1 lw 2 notitle, \
"../work_area_poisson_test_G2/plot_averaged_error_mpt1.dat" \
u ( 1.25+($1)):3 w l lt 2 lw 2 notitle, \
"../work_area_poisson_test_G2/plot_averaged_error_mpt2.dat" \
u (-1.25+($1)):3 w l lt 2 lw 2 notitle, \
"../work_area_poisson_test_G2/plot_averaged_error_mpt3.dat" \
u ( 0.00+($1)):3 w l lt 2 lw 2 notitle, \
"../work_area_poisson_test_G3/plot_averaged_error_mpt1.dat" \
u ( 1.25+($1)):3 w l lt 3 lw 2 notitle, \
"../work_area_poisson_test_G3/plot_averaged_error_mpt2.dat" \
u (-1.25+($1)):3 w l lt 3 lw 2 notitle, \
"../work_area_poisson_test_G3/plot_averaged_error_mpt3.dat" \
u ( 0.00+($1)):3 w l lt 3 lw 2 notitle, \
"../work_area_poisson_test_G4/plot_averaged_error_mpt1.dat" \
u ( 1.25+($1)):3 w l lt 4 lw 2 notitle, \
"../work_area_poisson_test_G4/plot_averaged_error_mpt2.dat" \
u (-1.25+($1)):3 w l lt 4 lw 2 notitle, \
"../work_area_poisson_test_G4/plot_averaged_error_mpt3.dat" \
u ( 0.00+($1)):3 w l lt 4 lw 2 notitle


#pause -1