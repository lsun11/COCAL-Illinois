set terminal postscript enhanced eps color dashed "Helvetica" 24
#set terminal postscript enhanced eps color solid "Helvetica" 16

set xlabel "{/Helvetica=38 x}"
set grid

set logscale x
set logscale y
set format x "10^{%L}"
set format y "10^{%L}"
set zero 1e-20
set mxtics 10 ; set mytics 10
set pointsize 1.6
set tics scale  2, 1

set key right top reverse Left
#set key width -1
set key box lw 0.2
set key title "df/dr 3rd"

set xrange [0.1:1000000]
set yrange [0.0001:10]

set ylabel "{/Helvetica=44 <}{/Helvetica=32 |}\
{/Symbol=32 dy}{/Helvetica=32 /}{/Symbol=32 y}\
{/Helvetica=32 |}{/Helvetica=44 > }{/Helvetica=32 [%]}"
set output 'plot_2pot_averaged_error_mpt_psi_logr.eps'
plot \
"../work_area_poisson_test_E1/plot_averaged_error_mpt1.dat" \
u ( 1.25+($1)):2 w l lt 1 lw 2 title "{E1}", \
"../work_area_poisson_test_E1/plot_averaged_error_mpt2.dat" \
u (-1.25+($1)):2 w l lt 1 lw 2 notitle, \
"../work_area_poisson_test_E2/plot_averaged_error_mpt1.dat" \
u ( 1.25+($1)):2 w l lt 2 lw 2 title "{E2}", \
"../work_area_poisson_test_E2/plot_averaged_error_mpt2.dat" \
u (-1.25+($1)):2 w l lt 2 lw 2 notitle, \
"../work_area_poisson_test_E3/plot_averaged_error_mpt1.dat" \
u ( 1.25+($1)):2 w l lt 3 lw 2 title "{E3}", \
"../work_area_poisson_test_E3/plot_averaged_error_mpt2.dat" \
u (-1.25+($1)):2 w l lt 3 lw 2 notitle, \
"../work_area_poisson_test_E4/plot_averaged_error_mpt1.dat" \
u ( 1.25+($1)):2 w l lt 4 lw 2 title "{E4}", \
"../work_area_poisson_test_E4/plot_averaged_error_mpt2.dat" \
u (-1.25+($1)):2 w l lt 4 lw 2 notitle


set label "COCP-1"at 0.13, 300 
set label "COCP-2"at 0.6, 50 

set yrange [0.0001:1000]
set ylabel "{/Helvetica=44 <}{/Helvetica=32 |}\
{/Symbol=32 da}{/Helvetica=32 /}{/Symbol=32 a}\
{/Helvetica=32 |}{/Helvetica=44 > }{/Helvetica=32 [%]}"
set output 'plot_2pot_averaged_error_mpt_alpha_logr.eps'
plot \
"../work_area_poisson_test_E1/plot_averaged_error_mpt1.dat" \
u ( 1.25+($1)):3 w l lt 1 lw 2 title "{E1}", \
"../work_area_poisson_test_E1/plot_averaged_error_mpt2.dat" \
u (-1.25+($1)):3 w l lt 1 lw 2 notitle, \
"../work_area_poisson_test_E2/plot_averaged_error_mpt1.dat" \
u ( 1.25+($1)):3 w l lt 2 lw 2 title "{E2}", \
"../work_area_poisson_test_E2/plot_averaged_error_mpt2.dat" \
u (-1.25+($1)):3 w l lt 2 lw 2 notitle, \
"../work_area_poisson_test_E3/plot_averaged_error_mpt1.dat" \
u ( 1.25+($1)):3 w l lt 3 lw 2 title "{E3}", \
"../work_area_poisson_test_E3/plot_averaged_error_mpt2.dat" \
u (-1.25+($1)):3 w l lt 3 lw 2 notitle, \
"../work_area_poisson_test_E4/plot_averaged_error_mpt1.dat" \
u ( 1.25+($1)):3 w l lt 4 lw 2 title "{E4}", \
"../work_area_poisson_test_E4/plot_averaged_error_mpt2.dat" \
u (-1.25+($1)):3 w l lt 4 lw 2 notitle


#pause -1
