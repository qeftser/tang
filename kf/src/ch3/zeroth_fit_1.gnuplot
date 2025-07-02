
set yrange [ -2 : 2 ]

plot "dat/zeroth_fit_1.dat" title "estimate" with lines,      \
     "dat/zeroth_poly_1.dat" title "measurement" with points, \
     x lc rgb "red" title "signal" with lines
