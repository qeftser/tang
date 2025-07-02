
#include "numerical_basics.h"
#include "recursive_least_squares_filtering.h"
#include <stdio.h>

int main(void) {

   r2lsfs state;
   init_recursive_second_order_least_squares_filter(&state);

   step_recursive_second_order_least_squares_filter(1.2,1,&state);
   printf("a0: %f\n",state.estimate_a0);
   printf("a1: %f\n",state.estimate_a1);
   printf("a2: %f\n",state.estimate_a2);

   step_recursive_second_order_least_squares_filter(0.2,2,&state);
   printf("a0: %f\n",state.estimate_a0);
   printf("a1: %f\n",state.estimate_a1);
   printf("a2: %f\n",state.estimate_a2);
   
   step_recursive_second_order_least_squares_filter(2.9,3,&state);
   printf("a0: %f\n",state.estimate_a0);
   printf("a1: %f\n",state.estimate_a1);
   printf("a2: %f\n",state.estimate_a2);

   step_recursive_second_order_least_squares_filter(2.1,4,&state);
   printf("a0: %f\n",state.estimate_a0);
   printf("a1: %f\n",state.estimate_a1);
   printf("a2: %f\n",state.estimate_a2);


   return 0;
}
