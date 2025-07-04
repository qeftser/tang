
set yrange [ 0 : 400000 ]

plot "dat/second_obj_m.dat" title "measurement" with points, \
     "dat/second_obj_p.dat" title "estimate" with lines, \
     400000 - 6000*x - 16.1*x*x title "real" with lines
