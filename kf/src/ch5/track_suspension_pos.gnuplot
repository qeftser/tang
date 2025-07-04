
set yrange [ -0.1 : 0.8 ]

plot "dat/suspension_kf.dat" using 2 with lines title "actual",       \
     "dat/suspension_kf.dat" using 3 with lines title "kf Q=0",       \
     "dat/suspension_kf.dat" using 4 with lines title "kf Q=0.00001", \
