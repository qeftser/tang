
set yrange [ -10 : 10 ]

set title "x velocity error over time"

plot "dat/ekf.dat" using 13 title "error" with lines, \
     "dat/ekf.dat" using 18 title "" with lines lc rgb "grey", \
     "dat/ekf.dat" using 19 title "" with lines lc rgb "grey"
