
set yrange [ 0.0 : 0.5 ]

plot "dat/second_kp22.dat" title "kalman gain p=INF" with lines, \
     "dat/second_kp22_p100.dat" title "kalman gain p=100" with lines, \
     "dat/second_kp22_p0.1.dat" title "kalman gain p=0.1" with lines, \
     "dat/second_kp22_p0.dat" title "kalman gain p=0.0" with points, \
     "dat/second_rp22.dat" title "least squares stddev gain" with lines
