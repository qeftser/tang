
set nokey
set title "monte-carlo error results with theoretical bounds"

set yrange [ -100 : 100 ]

plot "dat/second_order_err0.dat" with lines,                  \
     "dat/second_order_err1.dat" with lines,                  \
     "dat/second_order_err2.dat" with lines,                  \
     "dat/second_order_err3.dat" with lines,                  \
     "dat/second_order_err4.dat" with lines,                  \
      (3*(300*(x**2)-30*x+2)*50.0) / (10*x*(10*x+1)*(10*x+2)) with lines lc rgb "grey", \
     -(3*(300*(x**2)-30*x+2)*50.0) / (10*x*(10*x+1)*(10*x+2)) with lines lc rgb "grey"
      


