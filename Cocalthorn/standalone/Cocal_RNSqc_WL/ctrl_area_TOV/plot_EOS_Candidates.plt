#set terminal postscript eps color enhanced dashed "Helvetica" 24
#set output  "All_candidates.eps"

set key right bottom
set linestyle line 1 lt 2 lw 3
set key box

set format x "10^{%L}"
set format y "10^{%L}"
set ticscale 2
set mxtics 10
set mytics 10

set log x
set log y

set xlabel "{/Symbol=32 r} [g/cm^3]"
set ylabel "{/Helvetica=32 p({/Symbol=32 r})}"

set xrange [10.**13:2.*10.**15]
set yrange [10.**31:10.**36]

c=2.99792458*10.**10
f(x,k,gamma) = k*x**gamma
   
gamma_low = 1.35692395
gamma_high = 3.

#c1
k0_c1 = 3.99873692*10.**(-08)*c**2
k1_c1 = 2.238721139*10.**(-31)*c**2
rhodiv_c1 = 1.417289866*10.**14.
#c2b
k0_c2b = 3.99873692*10.**(-08)*c**2
k1_c2b = 1.12201845*10.**(-31)*c**2
rhodiv_c2b = 2.15795830*10.**14.
#c3
k0_c3 = 3.99873692*10.**(-08)*c**2
# wrong parameter k1_c3 = 2.81838293126443*10.**(-31)*c**2
k1_c3 = 7.079457843841396*10.**(-31)*c**2
rhodiv_c3 = 7.03317468*10.**13.
#c4
k0_c4 = 3.99873692*10.**(-08)*c**2
k1_c4 = 1.77827941004*10.**(-31)*c**2
rhodiv_c4 = 1.630497500126*10.**14.
#c5
k0_c5 = 3.99873692*10.**(-08)*c**2
# wrong parameter k1_c5 = 1.412537544623*10.**(-31)*c**2
k1_c5 = 2.81838293126443*10.**(-31)*c**2
rhodiv_c5 = 1.2319617560063*10.**14.

plot x < rhodiv_c1  ? f(x,k0_c1,gamma_low)  : f(x,k1_c1,gamma_high) \
     lw 4.0 lt  1 title "c1", \
     x < rhodiv_c2b ? f(x,k0_c2b,gamma_low) : f(x,k1_c2b,gamma_high) \
     lw 4.0 lt  2 title "c2b",\
     x < rhodiv_c3  ? f(x,k0_c3,gamma_low)  : f(x,k1_c3,gamma_high) \
     lw 4.0 lt  3 title "c3" ,\
     x < rhodiv_c4  ? f(x,k0_c4,gamma_low)  : f(x,k1_c4,gamma_high) \
     lw 4.0 lt  4 title "c4" ,\
     x < rhodiv_c5  ? f(x,k0_c5,gamma_low)  : f(x,k1_c5,gamma_high) \
     lw 4.0 lt  5 title "c5"

pause -1
