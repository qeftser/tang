
set title "y velocity estimate error"

plot "dat/ekf.dat" using 13 with lines title "error", \
     "dat/ekf.dat" using 20 with lines lc rgb "grey" title "", \
     "dat/ekf.dat" using 21 with lines lc rgb "grey" title ""
