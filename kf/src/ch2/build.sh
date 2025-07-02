#!/bin/sh

# you really should have used a makefile for that...

# >_<

echo "gcc -o fit_zero.e fit_zero.c"
gcc -o fit_zero.e fit_zero.c ../numerical_basics.c ../method_of_least_squares.c -lm
echo "gcc -o fit_one.e fit_one.c"
gcc -o fit_one.e fit_one.c ../numerical_basics.c ../method_of_least_squares.c -lm
echo "gcc -o fit_two.e fit_two.c"
gcc -o fit_two.e fit_two.c ../numerical_basics.c ../method_of_least_squares.c -lm
echo "mkdir dat"
mkdir -p dat
