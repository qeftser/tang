
set title "x velocity and estimate"

set yrange [ -16000 : -6100 ]

plot "dat/fading.dat" using 3 with lines title "velocity", \
     "dat/fading.dat" using 6 with lines title "estimate"
