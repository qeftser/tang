
set nokey
set title "monte-carlo error results with theoretical bounds"

set yrange [ -200 : 200 ]

plot "dat/second_order_err_v0.dat" with lines,                  \
     "dat/second_order_err_v1.dat" with lines,                  \
     "dat/second_order_err_v2.dat" with lines,                  \
     "dat/second_order_err_v3.dat" with lines,                  \
     "dat/second_order_err_v4.dat" with lines,                  \
      (12*(1600*(x**2)-300*x+11)*50.0) / (10*x*(100*(x**2)-1)*(100*(x**2)-4)*0.01) with lines lc rgb "grey", \
     -(12*(1600*(x**2)-300*x+11)*50.0) / (10*x*(100*(x**2)-1)*(100*(x**2)-4)*0.01) with lines lc rgb "grey"
      


