
set nokey
set title "monte-carlo error results with theoretical bounds"

set yrange [ -10 : 10 ]

plot "dat/first_order_err_v0.dat" with lines,                  \
     "dat/first_order_err_v1.dat" with lines,                  \
     "dat/first_order_err_v2.dat" with lines,                  \
     "dat/first_order_err_v3.dat" with lines,                  \
     "dat/first_order_err_v4.dat" with lines,                  \
      60 / (10*x * (100*x**2 - 1) * 0.01) with lines lc rgb "grey",  \
     -60 / (10*x * (100*x**2 - 1) * 0.01) with lines lc rgb "grey"
