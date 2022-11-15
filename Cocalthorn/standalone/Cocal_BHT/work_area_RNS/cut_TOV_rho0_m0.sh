# m0_1.8_TOV.txt  is the output of findstats.f90. Here we cut the columns that correspond to the 
# rest mass density and the rest mass and then we create a file with these two numbers at columns
# 38 and 20 so as to plot together with RNS output.

cut -c49-69 m0_1.8_TOV.txt  > a1

cut -c233-253 m0_1.8_TOV.txt > a2

echo 0.0 > a0

#     1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 
paste a0 a0 a0 a0 a0 a0 a0 a0 a0 a0 a0 a0 a0 a0 a0 a0 a0 a0 a0 a2 a0 a0 a0 a0 a0 a0 a0 a0 a0 a0 a0 a0 a0 a0 a0 a0 a0 a1 a0 a0 a0 a0 a0 a0 a0 a0 a0 a0  > m0_1.8_TOV_plot.txt
rm a0 a1 a2
