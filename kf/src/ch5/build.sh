#!/bin/sh

# you really should have used a makefile for that...

# >_<
# >_<
# >_<
# >_<

echo "gcc -o sin_data.e sin_data.c"
gcc -o sin_data.e sin_data.c ../numerical_basics.c ../polynomial_kalman_filters.c ../kalman_filters_in_a_nonpolynomial_world.c -lm
echo "gcc -o road_data.e road_data.c"
gcc -o road_data.e road_data.c ../numerical_basics.c ../polynomial_kalman_filters.c ../kalman_filters_in_a_nonpolynomial_world.c -lm
echo "mkdir dat"
mkdir -p dat
