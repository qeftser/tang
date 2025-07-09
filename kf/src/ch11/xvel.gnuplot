
set title "x velocity and estimate"

plot "dat/ekf.dat" using 4 with lines title "position", \
     "dat/ekf.dat" using 8 with lines title "estimate"
