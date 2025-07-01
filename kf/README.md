
# kf

Kalman filter section. Contains work produced from the examples and code given in the book
Fundamentals of Kalman Filtering: A Practical Approach Third Edition by Paul Zarchan and 
Howard Musoff. It is a pretty nice book :)

## numerical_basics.h

Contains code implimenting the needed functionality discussed in chapter one of the book.
This contains a simple matrix library with the nessesary routines, as well as several
random number generators to aid in the generation of gaussian and white noise. The
matrix library is implimented kind of interestingly, but works well if you know
how to use it... :O

## method_of_least_squares.h

Contains an implimentation of the least squares filter as described in chapter two. Given
an array of measurements and the corresponding time values they were taken at, a line
of best fit for any degree of polynomial desired can be computed - at least until the 
precision of the double type runs out :)

## /ch2

Replication of the examples presented in chapter two. Run build.sh to compile the data
generators, then run each program to produce it's data. The resulting data can be
visualized using gnuplot and the given scripts.

