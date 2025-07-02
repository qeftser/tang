
set yrange [ -2 : 2 ]

plot "dat/zero_order_fit0.dat" title "estimate" with lines,      \
     "dat/zero_order_poly0.dat" title "measurement" with points, \
     0 lc rgb "red" title "signal" with lines

