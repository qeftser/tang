
set yrange [ -2 : 18 ]

plot "dat/first_order_fit.dat" title "estimate" with lines,      \
     "dat/first_order_poly.dat" title "measurement" with points, \
     x + 3 lc rgb "red" title "signal" with lines

