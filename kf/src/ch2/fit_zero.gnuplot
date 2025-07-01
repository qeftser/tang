
plot "dat/fit_zero_measurements.dat" title "measurements",      \
     "dat/fit_zero_zeroth.dat" title "zero-order" with lines,   \
     "dat/fit_zero_first.dat" title "first-order" with lines,   \
     "dat/fit_zero_second.dat" title "second-order" with lines, \
     3 title "source" with lines dt 10 lc rgb "red"
