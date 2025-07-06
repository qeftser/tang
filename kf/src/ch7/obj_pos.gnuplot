
set yrange [ 0 : 200000 ]

plot "dat/object_falling.dat" using 2 with lines title "position", \
     "dat/object_tracking.dat" using 2 with lines title "estimate"
