
plot "dat/fit_one_measurements.dat" title "measurement"        \
     "dat/fit_one_zeroth.dat" title "zero-order" with lines,   \
     "dat/fit_one_first.dat" title "first-order" with lines,   \
     "dat/fit_one_second.dat" title "second-order" with lines, \
     x + 3 title "source" with lines dt 10 lc rgb "red"
