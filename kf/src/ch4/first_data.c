
/* generate the data for examining the behavior of the 
 * first order polynomial kalman filter. What is interesting
 * here is that we don't actually care about the data that
 * is being fed in, just the behavior of the gain matrices
 */

#include "../numerical_basics.h"
#include "../polynomial_kalman_filters.h"
#include <stdio.h>

double p11(double k, double sigma) {
   return (2*(2*k-1)*sigma*sigma) / (k*(k+1));
}

double p22(double k, double sigma) {
   return (12*sigma*sigma) / (k*(k*k-1)*1);
}

int main(void) {

   SRAND();

   FILE * f_0 = fopen("dat/first_kp11.dat","w");
   FILE * f_1 = fopen("dat/first_kp22.dat","w");
   FILE * d_0 = fopen("dat/first_rp11.dat","w");
   FILE * d_1 = fopen("dat/first_rp22.dat","w");

   fopkfs state0;
   init_first_order_polynomial_kalman_filter(&state0);
   state0.R = 1.0;
   state0.Q = 0.0;
   state0.H->vals[0] = 1.0;
   state0.Phi->vals[0] = 1;
   state0.Phi->vals[1] = 1;
   state0.Phi->vals[3] = 1;
   state0.P->vals[0] = 99999999999999;
   state0.P->vals[3] = 99999999999999;


   for (int i = 1; i < 100; ++i) {

      step_first_order_polynomial_kalman_filter(1.0,1.0,&state0);
      fprintf(f_0,"%lf  %lf\n",(double)i,state0.P->vals[0]);
      fprintf(d_0,"%lf  %lf\n",(double)i,p11(i,1.0));
      fprintf(f_1,"%lf  %lf\n",(double)i,state0.P->vals[3]);
      fprintf(d_1,"%lf  %lf\n",(double)i,p22(i,1.0));

   }

   fclose(f_0);
   fclose(f_1);
   fclose(d_0);
   fclose(d_1);

   destroy_first_order_polynomial_kalman_filter(&state0);

   fopkfs state1;
   fopkfs state2;
   fopkfs state3;
   init_first_order_polynomial_kalman_filter(&state1);
   init_first_order_polynomial_kalman_filter(&state2);
   init_first_order_polynomial_kalman_filter(&state3);
   state1.R = state2.R = state3.R = 1.0;
   state1.Q = state2.Q = state3.Q = 0.0;
   state1.H->vals[0] = state2.H->vals[0] = state3.H->vals[0] = 1.0;
   state1.Phi->vals[0] = state2.Phi->vals[0] = state3.Phi->vals[0] = 1;
   state1.Phi->vals[1] = state2.Phi->vals[1] = state3.Phi->vals[1] = 1;
   state1.Phi->vals[3] = state2.Phi->vals[3] = state3.Phi->vals[3] = 1;
   state1.P->vals[0] = 100;
   state2.P->vals[0] = 0.1;
   state3.P->vals[0] = 0;
   state1.P->vals[3] = 100;
   state2.P->vals[3] = 0.1;
   state3.P->vals[3] = 0;

   FILE * f_2 = fopen("dat/first_kp11_p100.dat","w");
   FILE * f_3 = fopen("dat/first_kp11_p0.1.dat","w");
   FILE * f_4 = fopen("dat/first_kp11_p0.dat","w");
   FILE * f_5 = fopen("dat/first_kp22_p100.dat","w");
   FILE * f_6 = fopen("dat/first_kp22_p0.1.dat","w");
   FILE * f_7 = fopen("dat/first_kp22_p0.dat","w");

   for (int i = 1; i < 100; ++i) {

      step_first_order_polynomial_kalman_filter(1.0,1.0,&state1);
      step_first_order_polynomial_kalman_filter(1.0,1.0,&state2);
      step_first_order_polynomial_kalman_filter(1.0,1.0,&state3);
      fprintf(f_2,"%lf  %lf\n",(double)i,state1.P->vals[0]);
      fprintf(f_5,"%lf  %lf\n",(double)i,state1.P->vals[3]);
      fprintf(f_3,"%lf  %lf\n",(double)i,state2.P->vals[0]);
      fprintf(f_6,"%lf  %lf\n",(double)i,state2.P->vals[3]);
      fprintf(f_4,"%lf  %lf\n",(double)i,state3.P->vals[0]);
      fprintf(f_7,"%lf  %lf\n",(double)i,state3.P->vals[3]);
   }

   fclose(f_2);
   fclose(f_3);
   fclose(f_4);
   fclose(f_5);
   fclose(f_6);
   fclose(f_7);

   destroy_first_order_polynomial_kalman_filter(&state1);
   destroy_first_order_polynomial_kalman_filter(&state2);
   destroy_first_order_polynomial_kalman_filter(&state3);

   /* we are going to try and estimate the falling object
    * here - expect the filter to fail! */

   fopkfs state4;
   fopkfs state5;
   fopkfs state6;
   init_first_order_polynomial_kalman_filter(&state4);
   init_first_order_polynomial_kalman_filter(&state5);
   init_first_order_polynomial_kalman_filter(&state6);

   FILE * c_d = fopen("dat/first_obj_m.dat","w");
   FILE * c_0 = fopen("dat/first_obj_p.dat","w");
   FILE * c_1 = fopen("dat/first_obj_v.dat","w");
   FILE * c_2 = fopen("dat/first_obj_pn.dat","w");
   FILE * c_3 = fopen("dat/first_obj_vn.dat","w");
   FILE * c_4 = fopen("dat/first_obj_pg.dat","w");
   FILE * c_5 = fopen("dat/first_obj_vg.dat","w");

   double ts = 0.1;
   double sigma = 1000.0;
   state4.P->vals[0] = state5.P->vals[0] = state6.P->vals[0] = 999999999999999;
   state4.P->vals[3] = state5.P->vals[3] = state6.P->vals[3] = 999999999999999;
   state4.Phi->vals[0] = state5.Phi->vals[0] = state6.Phi->vals[0] = 1.0;
   state4.Phi->vals[1] = state5.Phi->vals[1] = state6.Phi->vals[1] =  ts;
   state4.Phi->vals[3] = state5.Phi->vals[3] = state6.Phi->vals[3] = 1.0;
   state4.H->vals[0] = state5.H->vals[0] = state6.H->vals[0] = 1.0;
   state4.R = state5.R = state6.R = sigma*sigma;
   state4.Q = state6.Q = 0.0;
   state5.Q = 10000;

   for (int i = 1; i <= 300; ++i) {

      double t = i*ts;
      double x = 400000 - 6000*t - (32.2*t*t)/2.0;
      double xd = -6000 - 32.2*t;
      double xn = x + grand(0.0,sigma);
      step_first_order_polynomial_kalman_filter(xn,ts,&state4);
      step_first_order_polynomial_kalman_filter(xn,ts,&state5);
      step_gravity_informed_first_order_polynomial_kalman_filter(xn,ts,&state6);

      fprintf(c_d,"%lf  %lf\n",t,xn);
      fprintf(c_0,"%lf  %lf\n",t,state4.A->vals[0] - x);
      fprintf(c_1,"%lf  %lf\n",t,state4.A->vals[1] - xd);
      fprintf(c_2,"%lf  %lf\n",t,state5.A->vals[0] - x);
      fprintf(c_3,"%lf  %lf\n",t,state5.A->vals[1] - xd);
      fprintf(c_4,"%lf  %lf\n",t,state6.A->vals[0] - x);
      fprintf(c_5,"%lf  %lf\n",t,state6.A->vals[1] - xd);

   }

   fclose(c_d);
   fclose(c_0);
   fclose(c_1);
   fclose(c_2);
   fclose(c_3);
   fclose(c_4);
   fclose(c_5);
   
   destroy_first_order_polynomial_kalman_filter(&state4);
   destroy_first_order_polynomial_kalman_filter(&state5);
   destroy_first_order_polynomial_kalman_filter(&state6);

   return 0;
}
