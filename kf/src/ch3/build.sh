#!/bin/sh

# you really should have used a makefile for that...

# >_<
# >_<

echo "gcc -o zeroth_data.e zeroth_data.c"
gcc -o zeroth_data.e zeroth_data.c ../numerical_basics.c ../recursive_least_squares_filtering.c -lm
echo "gcc -o first_data.e first_data.c"
gcc -o first_data.e first_data.c ../numerical_basics.c ../recursive_least_squares_filtering.c -lm
echo "gcc -o second_data.e second_data.c"
gcc -o second_data.e second_data.c ../numerical_basics.c ../recursive_least_squares_filtering.c -lm
echo "mkdir dat"
mkdir -p dat
