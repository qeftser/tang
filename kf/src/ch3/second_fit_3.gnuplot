
plot "dat/second_order_fit_3.dat" title "estimate" with lines,      \
     4*(x**3) + 3*(x**2) + 2*x + 1 lc rgb "red" title "signal" with lines
