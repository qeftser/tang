
set yrange [ -200 : 200 ]

plot "dat/second_order_fit_vel.dat" title "estimate" with lines,      \
     10*x - 2 lc rgb "red" title "signal" with lines
