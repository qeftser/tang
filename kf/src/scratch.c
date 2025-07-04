
#include "numerical_basics.h"
#include "polynomial_kalman_filters.h"
#include "recursive_least_squares_filtering.h"
#include <stdio.h>
#include <float.h>

int main(void) {

   double stddev = 50.0;
   double _infinity = 1e200;

   sopkfs state;
   init_second_order_polynomial_kalman_filter(&state);
   state.Phi->vals[0] = 1.0;
   state.Phi->vals[1] = 1.0;
   state.Phi->vals[2] = 0.5;
   state.Phi->vals[4] = 1.0;
   state.Phi->vals[5] = 1.0;
   state.Phi->vals[8] = 1.0;
   state.H->vals[0] = 1.0;
   state.R = stddev;
   state.Q = 0.0;
   state.P->vals[0] = _infinity;
   state.P->vals[4] = _infinity;
   state.P->vals[8] = _infinity;

   r2lsfs state2;
   init_recursive_second_order_least_squares_filter(&state2);

   step_second_order_polynomial_kalman_filter(1.2,1,&state);
   step_recursive_second_order_least_squares_filter(1.2,1,&state2);
   printf("a0: %f\n",state.A->vals[0]);
   printf("b0: %f\n",state2.estimate_a0);
   printf("a1: %f\n",state.A->vals[1]);
   printf("b1: %f\n",state2.estimate_a1);
   printf("a2: %f\n",state.A->vals[2]);
   printf("b2: %f\n",state2.estimate_a2);

   step_second_order_polynomial_kalman_filter(0.2,1,&state);
   step_recursive_second_order_least_squares_filter(0.2,2,&state2);
   printf("a0: %f\n",state.A->vals[0]);
   printf("b0: %f\n",state2.estimate_a0);
   printf("a1: %f\n",state.A->vals[1]);
   printf("b1: %f\n",state2.estimate_a1);
   printf("a2: %f\n",state.A->vals[2]);
   printf("b2: %f\n",state2.estimate_a2);
   
   step_second_order_polynomial_kalman_filter(2.9,1,&state);
   step_recursive_second_order_least_squares_filter(2.9,3,&state2);
   printf("a0: %f\n",state.A->vals[0]);
   printf("b0: %f\n",state2.estimate_a0);
   printf("a1: %f\n",state.A->vals[1]);
   printf("b1: %f\n",state2.estimate_a1);
   printf("a2: %f\n",state.A->vals[2]);
   printf("b2: %f\n",state2.estimate_a2);

   step_second_order_polynomial_kalman_filter(2.1,1,&state);
   step_recursive_second_order_least_squares_filter(2.1,4,&state2);
   printf("a0: %f\n",state.A->vals[0]);
   printf("b0: %f\n",state2.estimate_a0);
   printf("a1: %f\n",state.A->vals[1]);
   printf("b1: %f\n",state2.estimate_a1);
   printf("a2: %f\n",state.A->vals[2]);
   printf("b2: %f\n",state2.estimate_a2);

   destroy_second_order_polynomial_kalman_filter(&state);
   return 0;
}
