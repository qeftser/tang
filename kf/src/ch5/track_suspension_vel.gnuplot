

set yrange [ -2 : 4 ]

plot "dat/suspension_kf.dat" using 5 with lines title "actual",       \
     "dat/suspension_kf.dat" using 6 with lines title "kf Q=0",       \
     "dat/suspension_kf.dat" using 7 with lines title "kf Q=0.00001", \
