#!/bin/sh

# you really should have used a makefile for that...

# >_<
# >_<
# >_<
# >_<
# >_<
# >_<
# >_<
# >_<
# >_<

echo "gcc -o fading.e fading.c"
gcc -g -o fading.e fading.c ../numerical_basics.c -lm
echo "mkdir dat"
mkdir -p dat
