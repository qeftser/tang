
plot "dat/sin_track_pos_eq.dat" with lines title "first_order laplace", \
     "dat/sin_track_pos_t2.dat" with lines title "first_order taylor x2", \
     sin(x) with lines title "signal"
