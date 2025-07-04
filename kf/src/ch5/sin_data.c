
/* produce all the data related to the attempts to track a sin
 * wave using a kalman filter. */

#include "../kalman_filters_in_a_nonpolynomial_world.h"
#include "../polynomial_kalman_filters.h"
#include "../numerical_basics.h"
#include <stdio.h>
#include <math.h>

int main(void) {

   SRAND();

   double ts = 0.1;
   double _infinity = 999999999999999;

   fopkfs s0;
   fopkfs s1;
   sopkfs s2;
   sopkfs s3;
   fopkfs s4;
   fopkfs s5;

   init_first_order_polynomial_kalman_filter(&s0);
   init_first_order_polynomial_kalman_filter(&s1);
   init_second_order_polynomial_kalman_filter(&s2);
   init_second_order_polynomial_kalman_filter(&s3);
   init_first_order_polynomial_kalman_filter(&s4);
   init_first_order_polynomial_kalman_filter(&s5);

   s0.Phi->vals[0] = 1.0;
   s0.Phi->vals[1] = ts;
   s0.Phi->vals[3] = 1.0;
   s0.H->vals[0] = 1.0;
   s0.Q = 0.0;
   s0.R = 1.0;
   s0.P->vals[0] = _infinity;
   s0.P->vals[3] = _infinity;

   s1.Phi->vals[0] = 1.0;
   s1.Phi->vals[1] = ts;
   s1.Phi->vals[3] = 1.0;
   s1.H->vals[0] = 1.0;
   s1.Q = 10.0;
   s1.R = 1.0;
   s1.P->vals[0] = _infinity;
   s1.P->vals[3] = _infinity;

   s2.Phi->vals[0] = 1.0;
   s2.Phi->vals[1] = ts;
   s2.Phi->vals[2] = 0.5*ts*ts;
   s2.Phi->vals[4] = 1.0;
   s2.Phi->vals[5] = ts;
   s2.Phi->vals[8] = 1.0;
   s2.H->vals[0] = 1.0;
   s2.Q = 0.0;
   s2.R = 1.0;
   s2.P->vals[0] = _infinity;
   s2.P->vals[4] = _infinity;
   s2.P->vals[8] = _infinity;

   s3.Phi->vals[0] = 1.0;
   s3.Phi->vals[1] = ts;
   s3.Phi->vals[2] = 0.5*ts*ts;
   s3.Phi->vals[4] = 1.0;
   s3.Phi->vals[5] = ts;
   s3.Phi->vals[8] = 1.0;
   s3.H->vals[0] = 1.0;
   s3.Q = 10.0;
   s3.R = 1.0;
   s3.P->vals[0] = _infinity;
   s3.P->vals[4] = _infinity;
   s3.P->vals[8] = _infinity;

   s4.Phi->vals[0] = cos(ts);
   s4.Phi->vals[1] = sin(ts);
   s4.Phi->vals[2] = -sin(ts);
   s4.Phi->vals[3] = cos(ts);
   s4.H->vals[0] = 1.0;
   s4.Q = 0.0;
   s4.R = 1.0;
   s4.P->vals[0] = _infinity;
   s4.P->vals[3] = _infinity;

   s5.Phi->vals[0] = 1;
   s5.Phi->vals[1] = ts;
   s5.Phi->vals[2] = -ts;
   s5.Phi->vals[3] = 1;
   s5.H->vals[0] = 1.0;
   s5.Q = 0.0;
   s5.R = 1.0;
   s5.P->vals[0] = _infinity;
   s5.P->vals[3] = _infinity;

   FILE * f_d = fopen("dat/sin_data.dat","w");
   FILE * t_0 = fopen("dat/sin_track_pos.dat","w");
   FILE * t_1 = fopen("dat/sin_track_vel.dat","w");
   FILE * t_2 = fopen("dat/sin_track_pos_10.dat","w");
   FILE * t_3 = fopen("dat/sin_track_vel_10.dat","w");
   FILE * t_4 = fopen("dat/sin_track_pos_2nd.dat","w");
   FILE * t_5 = fopen("dat/sin_track_vel_2nd.dat","w");
   FILE * t_6 = fopen("dat/sin_track_pos_2nd_10.dat","w");
   FILE * t_7 = fopen("dat/sin_track_vel_2nd_10.dat","w");
   FILE * t_8 = fopen("dat/sin_track_pos_eq.dat","w");
   FILE * t_9 = fopen("dat/sin_track_vel_eq.dat","w");
   FILE * t_a = fopen("dat/sin_track_pos_t2.dat","w");
   FILE * t_b = fopen("dat/sin_track_vel_t2.dat","w");

   /* I ran the simulation for a lot longer because I was 
    * curious how the signal estimates would change over
    * time. The normal ones drop pretty quickly, but adding
    * noise seems to change the kalman gain equation quite
    * a lot. I am not sure what to think about that.
    */
   for (int i = 1; i <= 600; ++i) {

      double t = i*ts;
      double x = sin(t);
      double xn = x + grand(0.0,1.0);

      step_first_order_polynomial_kalman_filter(xn,ts,&s0);
      step_first_order_polynomial_kalman_filter(xn,ts,&s1);
      step_second_order_polynomial_kalman_filter(xn,ts,&s2);
      step_second_order_polynomial_kalman_filter(xn,ts,&s3);
      step_first_order_polynomial_kalman_filter(xn,ts,&s4);
      step_first_order_polynomial_kalman_filter(xn,ts,&s5);

      fprintf(f_d,"%lf  %lf\n",t,xn);
      fprintf(t_0,"%lf  %lf\n",t,s0.A->vals[0]);
      fprintf(t_1,"%lf  %lf\n",t,s0.A->vals[1]);
      fprintf(t_2,"%lf  %lf\n",t,s1.A->vals[0]);
      fprintf(t_3,"%lf  %lf\n",t,s1.A->vals[1]);
      fprintf(t_4,"%lf  %lf\n",t,s2.A->vals[0]);
      fprintf(t_5,"%lf  %lf\n",t,s2.A->vals[1]);
      fprintf(t_6,"%lf  %lf\n",t,s3.A->vals[0]);
      fprintf(t_7,"%lf  %lf\n",t,s3.A->vals[1]);
      fprintf(t_8,"%lf  %lf\n",t,s4.A->vals[0]);
      fprintf(t_9,"%lf  %lf\n",t,s4.A->vals[1]);
      fprintf(t_a,"%lf  %lf\n",t,s5.A->vals[0]);
      fprintf(t_b,"%lf  %lf\n",t,s5.A->vals[1]);

   }

   fclose(f_d);
   fclose(t_0);
   fclose(t_1);
   fclose(t_2);
   fclose(t_3);
   fclose(t_4);
   fclose(t_5);
   fclose(t_6);
   fclose(t_7);
   fclose(t_8);
   fclose(t_9);
   fclose(t_a);
   fclose(t_b);

   destroy_first_order_polynomial_kalman_filter(&s0);
   destroy_first_order_polynomial_kalman_filter(&s1);
   destroy_second_order_polynomial_kalman_filter(&s2);
   destroy_second_order_polynomial_kalman_filter(&s3);
   destroy_first_order_polynomial_kalman_filter(&s4);
   destroy_first_order_polynomial_kalman_filter(&s5);

   return 0;
}
