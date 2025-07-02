
#include "../recursive_least_squares_filtering.h"
#include "../numerical_basics.h"
#include <stdlib.h>
#include <stdio.h>

int main(void) {

   SRAND();

   int monte_carlo_runs = 5;

   char * f_names[5] = { "dat/first_order_err0.dat",
                         "dat/first_order_err1.dat",
                         "dat/first_order_err2.dat",
                         "dat/first_order_err3.dat",
                         "dat/first_order_err4.dat"  };

   char * v_names[5] = { "dat/first_order_err_v0.dat",
                         "dat/first_order_err_v1.dat",
                         "dat/first_order_err_v2.dat",
                         "dat/first_order_err_v3.dat",
                         "dat/first_order_err_v4.dat"  };

   int num_samples = 100;
   double timestep = 0.1;

   for (int i = 0; i < monte_carlo_runs; ++i) {

      /* collect and record the measurement data for the zero order 
       * polynomial to fit.                                        */
      double * x_0 = malloc(sizeof(double) * num_samples);
      double * y_0 = malloc(sizeof(double) * num_samples);

      for (int i = 1; i <= num_samples; ++i) {
         x_0[i] = timestep * i;
         y_0[i] = x_0[i] + 3 + 5.0 * grand(0.0,1.0);
      }

      /* perform recursive least squares filtering on the data, 
       * writing each step to a file for plotting. */
      FILE * f_0 = fopen(f_names[i],"w+");
      FILE * v_0 = fopen(v_names[i],"w+");

      r1lsfs state;
      init_recursive_first_order_least_squares_filter(&state);

      for (int i = 0; i < num_samples; ++i) {

         step_recursive_first_order_least_squares_filter(y_0[i],x_0[i]+timestep,&state);

         double error = state.estimate_a0 - (x_0[i] + 3);

         fprintf(f_0,"%lf  %lf\n",x_0[i],error);
         fprintf(v_0,"%lf  %lf\n",x_0[i],1 - state.estimate_a1);

      }

      fclose(f_0);
      fclose(v_0);

      free(x_0);
      free(y_0);
   }


   /* now that we have our monte-carlo data, collect a single run of a 
    * noisy straight line to see the filter converge */

   FILE * f_1 = fopen("dat/first_order_fit.dat","w+");
   FILE * d_1 = fopen("dat/first_order_poly.dat","w+");

   double * x_1 = malloc(sizeof(double) * num_samples);
   double * y_1 = malloc(sizeof(double) * num_samples);

   r1lsfs state;
   init_recursive_first_order_least_squares_filter(&state);

   for (int i = 0; i < num_samples; ++i) {

      x_1[i] = i * timestep;
      y_1[i] = x_1[i] + 3 + 5 * grand(0.0,1.0);
      step_recursive_first_order_least_squares_filter(y_1[i],x_1[i]+timestep,&state);

      fprintf(d_1,"%lf  %lf\n",x_1[i],y_1[i]);
      fprintf(f_1,"%lf  %lf\n",x_1[i],state.estimate_a0 + state.estimate_a1 * 0.1);

   }

   fclose(d_1);
   fclose(f_1);
   free(x_1);
   free(y_1);

   /* now pit the filter against a noiseless quadratic to watch it diverge */

   FILE * f_2 = fopen("dat/first_order_fit_2.dat","w+");

   double * x_2 = malloc(sizeof(double) * num_samples);
   double * y_2 = malloc(sizeof(double) * num_samples);

   init_recursive_first_order_least_squares_filter(&state);

   for (int i = 0; i < num_samples; ++i) {

      x_2[i] = i * timestep;
      y_2[i] = 3 * x_2[i] * x_2[i] + 2 * x_2[i] + 1;
      step_recursive_first_order_least_squares_filter(y_2[i],x_2[i]+timestep,&state);

      fprintf(f_2,"%lf  %lf\n",x_2[i],state.estimate_a0);

   }

   fclose(f_2);
   free(x_2);
   free(y_2);

   return 0;
}
