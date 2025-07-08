
set title "y position error over time"

plot "dat/ekf.dat" using 14 title "error" with lines, \
     "dat/ekf.dat" using 20 title "" with lines lc rgb "grey", \
     "dat/ekf.dat" using 21 title "" with lines lc rgb "grey"
