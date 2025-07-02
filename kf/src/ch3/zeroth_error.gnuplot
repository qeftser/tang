
set nokey
set title "monte-carlo error results with theoretical bounds"

set yrange [ -2 : 2 ]

plot "dat/zero_order_fit0.dat" with lines,  \
     "dat/zero_order_fit1.dat" with lines,  \
     "dat/zero_order_fit2.dat" with lines,  \
     "dat/zero_order_fit3.dat" with lines,  \
     "dat/zero_order_fit4.dat" with lines,  \
     0 lc rgb "red",                        \
      1 / x with lines lc rgb "grey",       \
     -1 / x with lines lc rgb "grey"

