
set yrange [ -2 : 2 ]

plot "dat/sin_track_vel_eq.dat" with lines title "first_order laplace", \
     "dat/sin_track_vel_t2.dat" with lines title "first_order taylor x2", \
     cos(x) with lines title "signal"
