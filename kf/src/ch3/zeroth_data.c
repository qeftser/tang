
#include "../recursive_least_squares_filtering.h"
#include "../numerical_basics.h"
#include <stdlib.h>
#include <stdio.h>

int main(void) {

   SRAND();

   int monte_carlo_runs = 5;

   char * d_names[5] = { "dat/zero_order_poly0.dat",
                         "dat/zero_order_poly1.dat",
                         "dat/zero_order_poly2.dat",
                         "dat/zero_order_poly3.dat",
                         "dat/zero_order_poly4.dat"  };

   char * f_names[5] = { "dat/zero_order_fit0.dat",
                         "dat/zero_order_fit1.dat",
                         "dat/zero_order_fit2.dat",
                         "dat/zero_order_fit3.dat",
                         "dat/zero_order_fit4.dat"  };

   int num_samples = 100;
   double timestep = 0.1;

   for (int i = 0; i < monte_carlo_runs; ++i) {

      /* collect and record the measurement data for the zero order 
       * polynomial to fit.                                        */
      double * x_0 = malloc(sizeof(double) * num_samples);
      double * y_0 = malloc(sizeof(double) * num_samples);

      for (int i = 0; i < num_samples; ++i) {
         x_0[i] = timestep * i;

         /* note: in the book they use 1 as the line to fit. Here I am using 
          *       0 because it keeps me from also having to generate error
          *       files for each run. The estimate runs will be the error :)
          */
         y_0[i] = 0 + grand(0.0,1.0);
      }

      FILE * d_0 = fopen(d_names[i],"w+");

      for (int i = 0; i < num_samples; ++i)
         fprintf(d_0,"%lf  %lf\n",x_0[i],y_0[i]);

      fclose(d_0);

      /* perform recursive least squares filtering on the data, 
       * writing each step to a file for plotting. */
      FILE * f_0 = fopen(f_names[i],"w+");

      r0lsfs state;
      init_recursive_zeroth_order_least_squares_filter(&state);

      for (int i = 0; i < num_samples; ++i) {

         step_recursive_zeroth_order_least_squares_filter(y_0[i],x_0[i],&state);

         fprintf(f_0,"%lf  %lf\n",x_0[i],state.estimate);

      }

      fclose(f_0);

      free(x_0);
      free(y_0);
   }


   /* now that we have our monte-carlo data, collect a single run of a 
    * noiseless straight line to watch the filter diverge. */

   FILE * d_1 = fopen("dat/zeroth_poly_1.dat","w+");
   FILE * f_1 = fopen("dat/zeroth_fit_1.dat","w+");

   double * x_1 = malloc(sizeof(double) * num_samples);
   double * y_1 = malloc(sizeof(double) * num_samples);

   r0lsfs state;
   init_recursive_zeroth_order_least_squares_filter(&state);

   for (int i = 0; i < num_samples; ++i) {

      x_1[i] = i * timestep;
      y_1[i] = x_1[i];
      step_recursive_zeroth_order_least_squares_filter(y_1[i],x_1[i],&state);

      fprintf(d_1,"%lf  %lf\n",x_1[i],y_1[i]);
      fprintf(f_1,"%lf  %lf\n",x_1[i],state.estimate);

   }

   fclose(d_1);
   fclose(f_1);
   free(x_1);
   free(y_1);


   return 0;
}
