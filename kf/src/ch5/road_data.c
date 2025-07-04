
/* replication of the results from the section on 
 * applying the kalman filter to a suspension system.
 *
 * The things my neighbor said about the engineer/compuer
 * programmer working with the mathmatician makes more
 * sense when I see this kind of stuff...
 */

#include "../numerical_basics.h"
#include "../polynomial_kalman_filters.h"
#include "../kalman_filters_in_a_nonpolynomial_world.h"
#include <stdio.h>
#include <math.h>

int main(void) {

   double w = 6.28, wn = 6.28*0.1;
   double t;
   double zeta = 0.7;
   double h = 0.001;
   double A = 0.1;
   double x1 = 0.25, x1d = 0;
   double x2, x2d, x2dd, x1dd, x1d_old, x1_old, x2dd_old = 0.0;
   double s = 0.0;
   double b1 = 1.0, b2 = 1.25;

   double _infinity = 999999999999999;
   double sigma = 0.1;

   SRAND();

   /* start with a program that performs the Runge-Kutta 
    * integration on the suspension system and displays
    * it's movements and the corresponding behavior of
    * the actual system
    */

   FILE * f_0 = fopen("dat/suspension.dat","w");
   FILE * f_1 = fopen("dat/height.dat","w");

   t = 0.0;
   while (t < 20) {
      s += h;

      x1_old = x1;
      x1d_old = x1d;

      x2 = A*sin(w*t);
      x2d = A*w*cos(w*t);
      x2dd = -A*w*w*sin(w*t);
      x1dd = -2*zeta*wn*x1d - wn*wn*x1 - x2dd;
      x1 = x1 + x1d*h;
      x1d = x1d + x1dd*h;

      t += h;

      x2 = A*sin(w*t);
      x2d = A*w*cos(w*t);
      x2dd = -A*w*w*sin(w*t);
      x1dd = -2*zeta*wn*x1d - wn*wn*x1 - x2dd;
      x1 = 0.5*(x1_old+x1+h*x1d);
      x1d = 0.5*(x1d_old+x1d+h*x1dd);

      if (s > 0.099999) {
         s = 0.0;
         fprintf(f_0,"%lf  %lf\n",t,x1 + b2);
         fprintf(f_1,"%lf  %lf\n",t,x1 + x2 + b1 + b2);
      }
   }

   fclose(f_0);
   fclose(f_1);

   /* now attempt to apply the kalman filter to the measurement
    * of the signal's position. This will be done with a filter
    * that has process noise as well as one that does not. 
    *
    * It seems that the theoretical performance bounds of the
    * kalman filter are actually just the diagonals of the process
    * noise matrix. This is nice, as they are easy to access and
    * print along with the data. This also makes sense, as they
    * are the expected covariance of the measurements, which is
    * just sigma^2 for the diagomals - which is just stddev :)
    */

   /* apparently you can store lots of data in a single file. I feel
    * kind of silly now.. */
   FILE * f_2 = fopen("dat/suspension_kf.dat","w");

   double a = -zeta*wn;
   double b = wn*sqrt(1 - zeta*zeta);

   ffopkfs state0;
   ffopkfs state1;
   init_full_first_order_polynomial_kalman_filter(&state0);
   init_full_first_order_polynomial_kalman_filter(&state1);

   state0.P->vals[0] = state1.P->vals[0] = _infinity;
   state0.P->vals[3] = state1.P->vals[3] = _infinity;
   state0.R = state1.R = sigma;
   state0.Q = 0.0;
   state1.Q = 0.00001;
   state0.H->vals[0] = state1.H->vals[0] = 1.0;
   state0.G->vals[0] = state1.G->vals[0] = -(exp(a*0.1)*(-a*sin(b*0.1)-b*cos(b*0.1)) + b)/(b*(a*a+b*b));
   state0.G->vals[1] = state1.G->vals[1] = -(exp(a*0.1)*sin(b*0.1))/b;
   state0.Phi->vals[0] = state1.Phi->vals[0] = (exp(a*0.1)*(-a*sin(b*0.1)+b*cos(b*0.1)))/b;
   state0.Phi->vals[1] = state1.Phi->vals[1] = (exp(a*0.1)*sin(b*0.1))/b;
   state0.Phi->vals[2] = state1.Phi->vals[2] = (-wn*wn*exp(a*0.1)*sin(b*0.1))/b;
   state0.Phi->vals[3] = state1.Phi->vals[3] = (exp(a*0.1)*(a*sin(b*0.1) + b*cos(b*0.1)))/b;

   t = 0.0;
   s = 0.0;
   x1 = 0.25; 
   x1d = 0;
   while (t < 20) {
      s += h;

      x1_old = x1;
      x1d_old = x1d;

      x2 = A*sin(w*t);
      x2d = A*w*cos(w*t);
      x2dd = -A*w*w*sin(w*t);
      x1dd = -2*zeta*wn*x1d - wn*wn*x1 - x2dd;
      x1 = x1 + x1d*h;
      x1d = x1d + x1dd*h;

      t += h;

      x2 = A*sin(w*t);
      x2d = A*w*cos(w*t);
      x2dd = -A*w*w*sin(w*t);
      x1dd = -2*zeta*wn*x1d - wn*wn*x1 - x2dd;
      x1 = 0.5*(x1_old+x1+h*x1d);
      x1d = 0.5*(x1d_old+x1d+h*x1dd);

      if (s >= 0.099999) {
         s = 0.0;

         double z = x1 + grand(0.0,sigma);

         step_full_first_order_polynomial_kalman_filter(z,x2dd_old,0.1,&state0);
         step_full_first_order_polynomial_kalman_filter(z,x2dd_old,0.1,&state1);

         fprintf(f_2,"%lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf"
                     "  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n",
                 t,x1,state0.A->vals[0],state1.A->vals[0],
                 x1d,state0.A->vals[1],state1.A->vals[1],
                 state0.A->vals[0] - x1, state1.A->vals[0] - x1,
                 state0.A->vals[1] - x1d, state1.A->vals[1] - x1d,
                 sqrt(state0.P->vals[0]),-sqrt(state0.P->vals[0]),
                 sqrt(state1.P->vals[0]),-sqrt(state1.P->vals[0]),
                 sqrt(state0.P->vals[1]),-sqrt(state0.P->vals[1]),
                 sqrt(state1.P->vals[1]),-sqrt(state1.P->vals[1]));

         x2dd_old = x2dd;
      }
   }

   destroy_full_first_order_polynomial_kalman_filter(&state0);
   destroy_full_first_order_polynomial_kalman_filter(&state1);

   fclose(f_2);

   return 0;
}
