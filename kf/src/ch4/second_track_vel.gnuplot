
set yrange [ -10000 : 0 ]

plot "dat/second_obj_v.dat" title "measurement" with lines, \
     -6000 - 32.2*x title "real" with lines
