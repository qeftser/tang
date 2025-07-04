#!/bin/sh

# you really should have used a makefile for that...

# >_<
# >_<
# >_<

echo "gcc -o first_data.e first_data.c"
gcc -o first_data.e first_data.c ../numerical_basics.c ../polynomial_kalman_filters.c -lm
echo "gcc -o second_data.e second_data.c"
gcc -o second_data.e second_data.c ../numerical_basics.c ../polynomial_kalman_filters.c -lm
#echo "gcc -o grav_data.e grav_data.c"
#echo "gcc -o accel_data.e accel_data.c"
echo "mkdir dat"
mkdir -p dat
