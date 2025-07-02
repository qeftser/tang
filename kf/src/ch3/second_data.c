
#include "../recursive_least_squares_filtering.h"
#include "../numerical_basics.h"
#include <stdlib.h>
#include <stdio.h>

int main(void) {

   SRAND();

   int monte_carlo_runs = 5;

   char * f_names[5] = { "dat/second_order_err0.dat",
                         "dat/second_order_err1.dat",
                         "dat/second_order_err2.dat",
                         "dat/second_order_err3.dat",
                         "dat/second_order_err4.dat"  };

   char * v_names[5] = { "dat/second_order_err_v0.dat",
                         "dat/second_order_err_v1.dat",
                         "dat/second_order_err_v2.dat",
                         "dat/second_order_err_v3.dat",
                         "dat/second_order_err_v4.dat"  };

   char * a_names[5] = { "dat/second_order_err_a0.dat",
                         "dat/second_order_err_a1.dat",
                         "dat/second_order_err_a2.dat",
                         "dat/second_order_err_a3.dat",
                         "dat/second_order_err_a4.dat"  };

   int num_samples = 100;
   double timestep = 0.1;

   for (int i = 0; i < monte_carlo_runs; ++i) {

      /* collect and record the measurement data for the zero order 
       * polynomial to fit.                                        */
      double * x_0 = malloc(sizeof(double) * num_samples);
      double * y_0 = malloc(sizeof(double) * num_samples);

      for (int i = 0; i < num_samples; ++i) {
         x_0[i] = timestep * i + timestep;
         y_0[i] = 5 * x_0[i] * x_0[i] - 2 * x_0[i] + 2 + 50.0 * grand(0.0,1.0);
      }

      /* perform recursive least squares filtering on the data, 
       * writing each step to a file for plotting. */
      FILE * f_0 = fopen(f_names[i],"w+");
      FILE * v_0 = fopen(v_names[i],"w+");
      FILE * a_0 = fopen(a_names[i],"w+");

      r2lsfs state;
      init_recursive_second_order_least_squares_filter(&state);

      for (int i = 0; i < num_samples; ++i) {

         step_recursive_second_order_least_squares_filter(y_0[i],x_0[i],&state);

         double error = state.estimate_a0 - (5 * x_0[i] * x_0[i] - 2 * x_0[i] + 2);

         fprintf(f_0,"%lf  %lf\n",x_0[i],error);
         fprintf(v_0,"%lf  %lf\n",x_0[i],state.estimate_a1 - (2 + 10 * x_0[i]));
         fprintf(a_0,"%lf  %lf\n",x_0[i],state.estimate_a2 - 10);

      }

      fclose(f_0);
      fclose(v_0);
      fclose(a_0);

      free(x_0);
      free(y_0);
   }


   /* now that we have our monte-carlo data, collect a single run of a 
    * noisy straight line to see the filter converge */

   FILE * f_1 = fopen("dat/second_order_fit.dat","w+");
   FILE * f_v = fopen("dat/second_order_fit_vel.dat","w+");
   FILE * f_a = fopen("dat/second_order_fit_acc.dat","w+");
   FILE * d_1 = fopen("dat/second_order_poly.dat","w+");

   double * x_1 = malloc(sizeof(double) * num_samples);
   double * y_1 = malloc(sizeof(double) * num_samples);

   r2lsfs state;
   init_recursive_second_order_least_squares_filter(&state);

   for (int i = 0; i < num_samples; ++i) {

      x_1[i] = i * timestep + timestep;
      y_1[i] = 5 * x_1[i] * x_1[i] - 2 * x_1[i] + 2 + 50.0 * grand(0.0,1.0);
      step_recursive_second_order_least_squares_filter(y_1[i],x_1[i],&state);

      fprintf(d_1,"%lf  %lf\n",x_1[i],y_1[i]);
      fprintf(f_1,"%lf  %lf\n",x_1[i],state.estimate_a0);
      fprintf(f_v,"%lf  %lf\n",x_1[i],state.estimate_a1);
      fprintf(f_a,"%lf  %lf\n",x_1[i],state.estimate_a2);

   }

   fclose(d_1);
   fclose(f_1);
   fclose(f_v);
   fclose(f_a);
   free(x_1);
   free(y_1);

   /* now pit the filter against a noiseless quadratic to watch it diverge */

   FILE * f_2 = fopen("dat/second_order_fit_3.dat","w+");

   double * x_2 = malloc(sizeof(double) * num_samples);
   double * y_2 = malloc(sizeof(double) * num_samples);

   init_recursive_second_order_least_squares_filter(&state);

   for (int i = 0; i < num_samples; ++i) {

      x_2[i] = i * timestep + timestep;
      y_2[i] = 4 * x_2[i] * x_2[i] * x_2[i] + 3 * x_2[i] * x_2[i] + 2 * x_2[i] + 1;
      step_recursive_second_order_least_squares_filter(y_2[i],x_2[i],&state);

      fprintf(f_2,"%lf  %lf\n",x_2[i],state.estimate_a0 + state.estimate_a1 * 0.1 + state.estimate_a2 * 0.01);

   }

   fclose(f_2);
   free(x_2);
   free(y_2);

   return 0;
}
