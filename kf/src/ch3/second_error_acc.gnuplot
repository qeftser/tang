
set nokey
set title "monte-carlo error results with theoretical bounds"

set yrange [ -50 : 50 ]

plot "dat/second_order_err_a0.dat" with lines,                  \
     "dat/second_order_err_a1.dat" with lines,                  \
     "dat/second_order_err_a2.dat" with lines,                  \
     "dat/second_order_err_a3.dat" with lines,                  \
     "dat/second_order_err_a4.dat" with lines,                  \
       (720*50.0) / (10*x*(100*(x**2)-1)*(100*(x**2)-4)*1.0e-4) with lines lc rgb "grey", \
      -(720*50.0) / (10*x*(100*(x**2)-1)*(100*(x**2)-4)*1.0e-4) with lines lc rgb "grey"
      


