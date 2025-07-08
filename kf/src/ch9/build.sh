#!/bin/sh

# you really should have used a makefile for that...

# >_<
# >_<
# >_<
# >_<
# >_<
# >_<
# >_<

echo "gcc -o ekf.e ekf.c"
gcc -g -o ekf.e ekf.c ../numerical_basics.c -lm
echo "gcc -o duel_kf.e duel_kf.c"
gcc -g -o duel_kf.e duel_kf.c ../numerical_basics.c -lm
echo "mkdir dat"
mkdir -p dat
