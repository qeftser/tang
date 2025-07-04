

set yrange [ -0.04 : 0.04 ]

plot "dat/suspension_kf.dat" using 10 with lines title "kf Q=0", \
     "dat/suspension_kf.dat" using 11 with lines title "kf Q=0.00001", \
     "dat/suspension_kf.dat" using 16 title "" with lines lc rgb "grey", \
     "dat/suspension_kf.dat" using 17 title "" with lines lc rgb "grey", \
     "dat/suspension_kf.dat" using 18 title "" with lines lc rgb "grey", \
     "dat/suspension_kf.dat" using 19 title "" with lines lc rgb "grey"
