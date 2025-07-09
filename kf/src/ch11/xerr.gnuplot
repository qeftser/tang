
set title "x position estimate error"

plot "dat/ekf.dat" using 10 with lines title "error", \
     "dat/ekf.dat" using 14 with lines lc rgb "grey" title "", \
     "dat/ekf.dat" using 15 with lines lc rgb "grey" title ""
