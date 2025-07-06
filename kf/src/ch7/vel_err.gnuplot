
set yrange [ -100 : 100 ]

plot "dat/object_tracking.dat" using 5 with lines title "error", \
     "dat/object_tracking.dat" using 8 with lines lc rgb "grey" title "", \
     "dat/object_tracking.dat" using 9 with lines lc rgb "grey" title ""
