
plot "dat/first_order_fit_2.dat" title "estimate" with lines,      \
     3 * x**2 + 2 * x + 1 lc rgb "red" title "signal" with lines
