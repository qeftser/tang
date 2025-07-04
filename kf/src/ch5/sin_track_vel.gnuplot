
set yrange [ -3 : 3 ]

plot "dat/sin_track_vel.dat" with lines title "first_order Q=0", \
     "dat/sin_track_vel_10.dat" with lines title "first_order Q=10", \
     "dat/sin_track_vel_2nd.dat" with lines title "second_order Q=0", \
     "dat/sin_track_vel_2nd_10.dat" with lines title "second_order Q=10", \
     cos(x) with lines title "signal"
