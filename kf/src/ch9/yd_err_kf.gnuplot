
set yrange [ -100 : 100 ]

set title "y velocity error over time"

plot "dat/duel_kf.dat" using 15 title "error" with lines, \
     "dat/duel_kf.dat" using 22 title "" with lines lc rgb "grey", \
     "dat/duel_kf.dat" using 23 title "" with lines lc rgb "grey"
