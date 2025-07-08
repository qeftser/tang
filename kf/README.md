
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

## recursive_least_squares_filtering.h

Contains implimentations of the zeroth, first, and second order recursive least
squares filters. These are pretty nice and basically do the job of the least
squares filter in real time.

## /ch3

Replication of the examples presented in chapter three. Run build.sh to compile
the data generators. Each of these will generate data used by the plotting
scripts. The fitting works, but the error margins on the error plots are either
wrong or the scaling of the filters are wrong somewhere. Either way, it is not
a big enough deal for me to care about right now :/. The filters work!

## polynomial_kalman_filters.h

Contains implimentations of the zeroth, first, and second order polynomial
kalman filters discussed in chapter 4. These are only partially compelete,
as they lack the deterministic input term G in the state propagation. These
filters are sufficient to replicate the results presented in chapter 4.

## /ch4

Replication of the examples in chapter four. Formula is the same as
before. Run build.sh and the generated commands to produce the data,
then use the plotting scripts to observe the resulting behaviors.

## kalman_filters_in_a_nonpolynomial_world.h

Contains an impliementation of a first order kalman filter that has
the G matrix and accepts a static u term. 

## /ch5 

Replication of the examples in chapter five. I will note that
the code for the suspension tracking is not written correctly,
and as a result the filter behaves rather poorly. This is
likely do to an error copying the matrix math over, but
I feel I understand the text and would perfer not to spend
several hours trying to debug the code. Say what you will!

## /ch7

Replication of the example presented in chapter seven. There is
no associated header file because the extended kalman filter is
programmed directly into the simulation. This is done because
several states needed to be recomputed at each interval, and I
would have used it once anyway. Plus it is good to have more
experience playing with the Riccati and update equations, right?

## /ch8

Attempted implimentation of the extended kalman filter covered in
chapter eight. I was able to get decent tracking with 1/beta as
a state and a process noise of 1, but still not as good as what 
was presented in the text. I went over the equations several times
but maybe missed something? Linear filter was not implimented because
I already spent a while on the nonlinear one :/

## /ch9

Implimentation fo the extended kalman fiter and duel linear kalman
filters from chapter 9. These things are tedious and annoying to
debug. Anyway, they are there for your viewing pleasure.
