
set yrange [ -0.02 : 0.02 ]

plot "dat/suspension_kf.dat" using 8 with lines title "kf Q=0", \
     "dat/suspension_kf.dat" using 9 with lines title "kf Q=0.00001", \
     "dat/suspension_kf.dat" using 12 title "" with lines lc rgb "grey", \
     "dat/suspension_kf.dat" using 13 title "" with lines lc rgb "grey", \
     "dat/suspension_kf.dat" using 14 title "" with lines lc rgb "grey", \
     "dat/suspension_kf.dat" using 15 title "" with lines lc rgb "grey"
