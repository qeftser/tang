
set nokey
set title "monte-carlo error results with theoretical bounds"

set yrange [ -10 : 10 ]

plot "dat/first_order_err0.dat" with lines,                  \
     "dat/first_order_err1.dat" with lines,                  \
     "dat/first_order_err2.dat" with lines,                  \
     "dat/first_order_err3.dat" with lines,                  \
     "dat/first_order_err4.dat" with lines,                  \
      10 * (20*x - 1) / (10*x * (10*x+1)) with lines lc rgb "grey", \
     -10 * (20*x - 1) / (10*x * (10*x+1)) with lines lc rgb "grey"

