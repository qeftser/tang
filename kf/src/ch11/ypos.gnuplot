
set title "y position and estimate"

set yrange [ -200 : 200 ]

plot "dat/ekf.dat" using 3 with lines title "position", \
     "dat/ekf.dat" using 7 with lines title "estimate"
