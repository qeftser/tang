#!/bin/sh

# you really should have used a makefile for that...

# >_<
# >_<
# >_<
# >_<
# >_<

echo "gcc -o drag_data.e drag_data.c"
gcc -g -o drag_data.e drag_data.c ../numerical_basics.c -lm
echo "mkdir dat"
mkdir -p dat
