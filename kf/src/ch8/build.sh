#!/bin/sh

# you really should have used a makefile for that...

# >_<
# >_<
# >_<
# >_<
# >_<
# >_<

echo "gcc -o obj_data.e obj_data.c"
gcc -g -o obj_data.e obj_data.c ../numerical_basics.c -lm
echo "mkdir dat"
mkdir -p dat
