
set title "x position, measurement, and estimate"

plot "dat/fading.dat" using 2 with lines title "position", \
     "dat/fading.dat" using 4 with lines title "measurement", \
     "dat/fading.dat" using 5 with lines title "estimate"
