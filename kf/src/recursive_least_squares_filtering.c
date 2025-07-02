
#include "recursive_least_squares_filtering.h"
#include <assert.h>
#include <float.h>

void init_recursive_zeroth_order_least_squares_filter(r0lsfs * state) {
   state->estimate = 0.0;
   state->num_samples  = 0;
}

void step_recursive_zeroth_order_least_squares_filter(double measurement, double elapsed, r0lsfs * state) {

   state->num_samples += 1;
   state->estimate = state->estimate + (1.0 / state->num_samples) * (measurement - state->estimate);

}

void init_recursive_first_order_least_squares_filter(r1lsfs * state) {
   state->estimate_a0 = 0.0;
   state->estimate_a1 = 0.0;
   state->elapsed     = 0.0;
   state->num_samples = 0;
}

void step_recursive_first_order_least_squares_filter(double measurement, double elapsed, r1lsfs * state) {

   assert(elapsed > state->elapsed && "given time must be greater than last sample time");

   state->num_samples += 1;

   double timestep = elapsed - state->elapsed;

   double k1 = 
      (2.0 * (2.0*state->num_samples - 1)) / (state->num_samples * (state->num_samples + 1));
   double k2 = 
      6.0 / (state->num_samples * (state->num_samples + 1) * timestep);
   double res = measurement - state->estimate_a0 - state->estimate_a1 * timestep;

   double new_a0 = state->estimate_a0 + state->estimate_a1 * timestep + k1 * res;
   double new_a1 = state->estimate_a1 + k2 * res;

   state->estimate_a0 = new_a0;
   state->estimate_a1 = new_a1;
   state->elapsed = elapsed;
}

void init_recursive_second_order_least_squares_filter(r2lsfs * state) {
   state->estimate_a0 = 0.0;
   state->estimate_a1 = 0.0;
   state->estimate_a2 = 0.0;
   state->elapsed     = 0.0;
   state->num_samples = 0;
}

void step_recursive_second_order_least_squares_filter(double measurement, double elapsed, r2lsfs * state) {

   assert(elapsed > state->elapsed && "given time must be greater than last sample time");

   int k = state->num_samples += 1;

   double timestep = elapsed - state->elapsed;

   double k1 = (3.0 * (3 * k*k - 3 * k + 2)) / (k * (k + 1) * (k + 2));
   double k2 = (18.0 * (2 * k - 1)) / (k * (k + 1) * (k + 2) * timestep);
   double k3 = (60.0 / (k * (k + 1) * (k + 2) * timestep*timestep));
   double res = measurement - state->estimate_a0 
                            - state->estimate_a1 * timestep 
                            - 0.5 * state->estimate_a2 * timestep * timestep;

   double new_a0 = state->estimate_a0 + state->estimate_a1 * timestep + 0.5 * state->estimate_a2 * timestep * timestep + k1 * res;
   double new_a1 = state->estimate_a1 + state->estimate_a2 * timestep + k2 * res;
   double new_a2 = state->estimate_a2 + k3 * res;

   state->estimate_a0 = new_a0;
   state->estimate_a1 = new_a1;
   state->estimate_a2 = new_a2;
   state->elapsed = elapsed;
}

