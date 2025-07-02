
set yrange [ -100 : 500 ]

plot "dat/second_order_fit.dat" title "estimate" with lines,      \
     "dat/second_order_poly.dat" title "measurement" with points, \
     5*(x**2) - 2*x + 2 lc rgb "red" title "signal" with lines
