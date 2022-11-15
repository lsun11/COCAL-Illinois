set terminal postscript enhanced eps color dashed "Helvetica" 24 \
fontfile "/usr/share/texmf/fonts/type1/bluesky/cm/cmsy10.pfb"

#set xlabel "{/Helvetica=28 x}"
#set ylabel "{/Symbol=28 y}"
set grid
set key left top reverse Left
set mxtics 5
set mytics 5
set xrange [0:1.2]
set output 'vep_iteration400_mpt1_f.eps'
plot "vep_iteration400_mpt1.txt" u 2:3 every :::0::0 w lp lw 2 lt 1 lc rgb "#0000FF" title "vepxf",\
     "vep_iteration400_mpt1.txt" u 2:4 every :::0::0 w lp lw 2 lt 1 lc rgb "#00FF00" title "vepyf",\
     "vep_iteration400_mpt1.txt" u 2:5 every :::0::0 w lp lw 2 lt 1 lc rgb "#FF0000" title "vepzf"

reset
set grid
set key left top reverse Left
set mxtics 5
set mytics 5
set xrange [0:1.2]
set output 'vep_iteration400_mpt1_g.eps'
plot "vep_iteration400_mpt1.txt" u 2:3 every :::1::1 w lp lw 2 lt 1 lc rgb "#00FF00" title "vepxg",\
     "vep_iteration400_mpt1.txt" u 2:4 every :::1::1 w lp lw 2 lt 1 lc rgb "#0000FF" title "vepyg",\
     "vep_iteration400_mpt1.txt" u 2:5 every :::1::1 w lp lw 2 lt 1 lc rgb "#FF0000" title "vepzg"

