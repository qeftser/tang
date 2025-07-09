
set title "x position and estimate"

set yrange [ 0 : 20000 ]

plot "dat/ekf.dat" using 2 with lines title "position", \
     "dat/ekf.dat" using 6 with lines title "estimate"

