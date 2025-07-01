
plot "dat/fit_two_measurements.dat" title "measurement",       \
     "dat/fit_two_zeroth.dat" title "zero-order" with lines,   \
     "dat/fit_two_first.dat" title "first-order" with lines,   \
     "dat/fit_two_second.dat" title "second-order" with lines, \
     5*x**2 - 2*x +2 title "source" with lines dt 10 lc rgb "red"
