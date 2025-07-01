
#include "../method_of_least_squares.h"
#include <stdio.h>

/* generate the data for the attempt to fit the zeroth, first, and second
 * order least square filters to the equation f(x) = x + 3
 *
 * data will be written to 
 * fit_one_measurements.dat
 * fit_one_zeroth.dat
 * fit_one_first.dat
 * fit_one_second.dat
 *
 * 100 measurements from 0 to 10 in intervals of 0.1 will be taken
 * gaussian noise of mean 0 and sigma 1 that is scaled by 5 will be
 * added to the function
 */

int main(void) {

   SRAND();

   const int num_samples = 100;
   const double sample_step = 0.1;

   double * x = malloc(sizeof(double) * num_samples);
   double * y = malloc(sizeof(double) * num_samples);
   double * a = malloc(sizeof(double) * 3);

   for (int i = 0; i < num_samples; ++i) {
      x[i] = sample_step * i;
      y[i] = x[i] + 3 + 5 * grand(0.0,1.0);
   }

   /* write out the measurements */
   FILE * measurements = fopen("dat/fit_one_measurements.dat","w+");

   for (int i = 0; i < num_samples; ++i) {
      fprintf(measurements,"%lf  %lf\n",x[i],y[i]);
   }

   fclose(measurements);

   least_squares_filter(x,y,num_samples,a,0);

   /* write out the least squares line of best fit 
    * for the zeroth order                         */
   FILE * zeroth = fopen("dat/fit_one_zeroth.dat","w+");

   for (int i = 0; i < num_samples; ++i) {
      fprintf(zeroth,"%lf  %lf\n",x[i],a[0]);
   }

   fclose(zeroth);

   least_squares_filter(x,y,num_samples,a,1);

   /* write out the least squares line of best fit
    * for the first order                          */
   FILE * first = fopen("dat/fit_one_first.dat","w+");

   for (int i = 0; i < num_samples; ++i) {
      fprintf(first,"%lf  %lf\n",x[i],a[0] + a[1] * x[i]);
   }

   fclose(first);

   least_squares_filter(x,y,num_samples,a,2);

   /* write out the least squares line of best fit
    * for the first order                          */
   FILE * second = fopen("dat/fit_one_second.dat","w+");

   for (int i = 0; i < num_samples; ++i) {
      fprintf(second,"%lf  %lf\n",x[i],a[0] + a[1] * x[i] + a[2] * x[i] * x[i]);
   }

   fclose(second);

   free(x);
   free(y);
   free(a);

   return 0;
}
