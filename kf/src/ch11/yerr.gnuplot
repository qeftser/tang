
set title "y position estimate error"

plot "dat/ekf.dat" using 11 with lines title "error", \
     "dat/ekf.dat" using 16 with lines lc rgb "grey" title "", \
     "dat/ekf.dat" using 17 with lines lc rgb "grey" title ""
