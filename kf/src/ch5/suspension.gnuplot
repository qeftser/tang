
set yrange [ -0.2 : 3.0 ]

plot "dat/suspension.dat" title "suspension" with lines, \
     "dat/height.dat" title "car height" with lines,     \
     0.1 * sin(6.28*x) title "road" with lines 
