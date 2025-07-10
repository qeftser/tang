
set yrange [ -500 : 500 ]

set title "error in x position estimate"

plot "dat/fading.dat" using 7 with lines title "error",          \
     "dat/fading.dat" using 9 with lines lc rgb "grey" title "", \
     "dat/fading.dat" using 10 with lines lc rgb "grey" title ""

