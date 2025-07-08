
set yrange [ -100 : 100 ]

set title "y velocity error over time"

plot "dat/ekf.dat" using 15 title "error" with lines, \
     "dat/ekf.dat" using 22 title "" with lines lc rgb "grey", \
     "dat/ekf.dat" using 23 title "" with lines lc rgb "grey"
