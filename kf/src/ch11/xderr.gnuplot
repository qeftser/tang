
set title "x velocity estimate error"

plot "dat/ekf.dat" using 12 with lines title "error", \
     "dat/ekf.dat" using 18 with lines lc rgb "grey" title "", \
     "dat/ekf.dat" using 19 with lines lc rgb "grey" title ""
