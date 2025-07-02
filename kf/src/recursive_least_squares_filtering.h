
/* the zeroth, first, and second order recursive least squares filters, as
 * presented in chapter three of the book. Each filter accepts a measurement
 * and the current time step at every interval, and updates it's internal
 * state with each measurement.
 */

#ifndef __TANG_RECURSIVE_LEAST_SQUARES_FILTERING__

#define __TANG_RECURSIVE_LEAST_SQUARES_FILTERING__
#include <stdint.h>

/* internal state of the recursive least squares filter of the 
 * zeroth order. The best estimate of the polynomial is given
 * by estimate, and num_samples tracks the number of samples
 * collected since the filter started collecting samples.
 */
typedef struct recursive_zeroth_order_least_squares_filter_state {
   double estimate;
   uint64_t num_samples;
} r0lsfs;

/* set the initial values of the recusive least squares filter of the
 * zeroth order.  */
void init_recursive_zeroth_order_least_squares_filter(r0lsfs * state);

/* collect a measurement and update the state of the recursive zeroth order
 * least squares filter. measurement corresponds to the measurement collected,
 * elapsed corresponds to the time the measurement was collected, and state
 * is the filter state.
 *
 * Note: elapsed must be a value greater than the elapsed recorded in state,
 * i.e. time must have passed between the last measurement and the current one.
 */
void step_recursive_zeroth_order_least_squares_filter(double measurement, double elapsed, r0lsfs * state);

/* internal state of the recursive least squares filter of the 
 * first order. The best estimate of the polynomial is given
 * by estimate_a0 and estimate_a1. num_samples tracks the number 
 * of samples collected since the filter started collecting samples.
 * elapsed holds the last sample time provided.
 */
typedef struct recursive_first_order_least_squares_filter_state {
   double estimate_a0;
   double estimate_a1;
   double elapsed;
   uint64_t num_samples;
} r1lsfs;

/* set the initial values of the recusive least squares filter of the
 * first order.  */
void init_recursive_first_order_least_squares_filter(r1lsfs * state);

/* collect a measurement and update the state of the recursive first order
 * least squares filter. measurement corresponds to the measurement collected,
 * elapsed corresponds to the time the measurement was collected, and state
 * is the filter state.
 *
 * Note: elapsed must be a value greater than the elapsed recorded in state,
 * i.e. time must have passed between the last measurement and the current one.
 */
void step_recursive_first_order_least_squares_filter(double measurement, double elapsed, r1lsfs * state);

/* internal state of the recursive least squares filter of the 
 * second order. The best estimate of the polynomial is given
 * by estimate_a0, estimate_a1, and estimate_a2. num_samples 
 * tracks the number of samples collected since the filter 
 * started collecting samples. elapsed holds the last sample
 * time provided.
 */
typedef struct recursive_second_order_least_squares_filter_state {
   double estimate_a0;
   double estimate_a1;
   double estimate_a2;
   double elapsed;
   uint64_t num_samples;
} r2lsfs;

/* set the initial values of the recusive least squares filter of the
 * second order.  */
void init_recursive_second_order_least_squares_filter(r2lsfs * state);

/* collect a measurement and update the state of the recursive second order
 * least squares filter. measurement corresponds to the measurement collected,
 * elapsed corresponds to the time the measurement was collected, and state
 * is the filter state.
 *
 * Note: elapsed must be a value greater than the elapsed recorded in state,
 * i.e. time must have passed between the last measurement and the current one.
 */
void step_recursive_second_order_least_squares_filter(double measurement, double elapsed, r2lsfs * state);

#endif
 
