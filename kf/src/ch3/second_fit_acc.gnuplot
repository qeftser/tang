
set yrange [ -50 : 50 ]

plot "dat/second_order_fit_acc.dat" title "estimate" with lines,      \
     10 lc rgb "red" title "signal" with lines
