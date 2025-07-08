

set title "x position error over time"

plot "dat/duel_kf.dat" using 12 title "error" with lines, \
     "dat/duel_kf.dat" using 16 title "" with lines lc rgb "grey", \
     "dat/duel_kf.dat" using 17 title "" with lines lc rgb "grey"
