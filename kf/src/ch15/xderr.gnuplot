
set yrange [ -30 : 30 ]

set title "error in x velocity estimate"

plot "dat/fading.dat" using 8 with lines title "error",          \
     "dat/fading.dat" using 11 with lines lc rgb "grey" title "", \
     "dat/fading.dat" using 12 with lines lc rgb "grey" title ""

