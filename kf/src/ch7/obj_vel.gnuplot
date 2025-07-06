
set yrange [ -7000 : -3000 ]

plot "dat/object_falling.dat" using 3 with lines title "velocity", \
     "dat/object_tracking.dat" using 4 with lines title "estimate"
