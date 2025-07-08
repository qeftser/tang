
set title "x position error over time"

plot "dat/ekf.dat" using 12 title "error" with lines, \
     "dat/ekf.dat" using 16 title "" with lines lc rgb "grey", \
     "dat/ekf.dat" using 17 title "" with lines lc rgb "grey"
