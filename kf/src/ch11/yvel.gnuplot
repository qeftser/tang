
set title "y velocity and estimate"

set yrange [ -30 : 30 ]

plot "dat/ekf.dat" using 5 with lines title "position", \
     "dat/ekf.dat" using 9 with lines title "estimate"
