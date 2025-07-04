
#include "../numerical_basics.h"
#include "../polynomial_kalman_filters.h"
#include <stdio.h>

double p11(double k, double sigma) {
   return (3*(3*k*k-3*k+2)*sigma*sigma)/(k*(k+1)*(k+2));
}

double p22(double k, double sigma) {
   return (12*(16*k*k-30*k+11)*sigma*sigma)/(k*(k*k-1)*(k*k-4));
}

double p33(double k, double sigma) {
   return (720*sigma*sigma)/(k*(k*k-1)*(k*k-4));
}

int main(void) {

   SRAND();

   FILE * f_0 = fopen("dat/second_kp11.dat","w");
   FILE * f_1 = fopen("dat/second_kp22.dat","w");
   FILE * f_2 = fopen("dat/second_kp33.dat","w");
   FILE * d_0 = fopen("dat/second_rp11.dat","w");
   FILE * d_1 = fopen("dat/second_rp22.dat","w");
   FILE * d_2 = fopen("dat/second_rp33.dat","w");

   sopkfs state0;
   init_second_order_polynomial_kalman_filter(&state0);
   state0.R = 1.0;
   state0.Q = 0.0;
   state0.H->vals[0] = 1.0;
   state0.Phi->vals[0] = 1;
   state0.Phi->vals[1] = 1;
   state0.Phi->vals[2] = 0.5;
   state0.Phi->vals[4] = 1;
   state0.Phi->vals[5] = 1;
   state0.Phi->vals[8] = 1;
   state0.P->vals[0] = 99999999999999;
   state0.P->vals[4] = 99999999999999;
   state0.P->vals[8] = 99999999999999;


   for (int i = 1; i < 100; ++i) {

      step_second_order_polynomial_kalman_filter(1.0,1.0,&state0);
      fprintf(f_0,"%lf  %lf\n",(double)i,state0.P->vals[0]);
      fprintf(d_0,"%lf  %lf\n",(double)i,p11(i,1.0));
      fprintf(f_1,"%lf  %lf\n",(double)i,state0.P->vals[4]);
      fprintf(d_1,"%lf  %lf\n",(double)i,p22(i,1.0));
      fprintf(f_2,"%lf  %lf\n",(double)i,state0.P->vals[8]);
      fprintf(d_2,"%lf  %lf\n",(double)i,p33(i,1.0));

   }

   fclose(f_0);
   fclose(f_1);
   fclose(f_2);
   fclose(d_0);
   fclose(d_1);
   fclose(d_2);

   destroy_second_order_polynomial_kalman_filter(&state0);

   sopkfs state1;
   sopkfs state2;
   sopkfs state3;
   init_second_order_polynomial_kalman_filter(&state1);
   init_second_order_polynomial_kalman_filter(&state2);
   init_second_order_polynomial_kalman_filter(&state3);
   state1.R = state2.R = state3.R = 1.0;
   state1.Q = state2.Q = state3.Q = 0.0;
   state1.H->vals[0] = state2.H->vals[0] = state3.H->vals[0] = 1.0;
   state1.Phi->vals[0] = state2.Phi->vals[0] = state3.Phi->vals[0] = 1;
   state1.Phi->vals[1] = state2.Phi->vals[1] = state3.Phi->vals[1] = 1;
   state1.Phi->vals[2] = state2.Phi->vals[2] = state3.Phi->vals[2] = 0.5;
   state1.Phi->vals[4] = state2.Phi->vals[4] = state3.Phi->vals[4] = 1;
   state1.Phi->vals[5] = state2.Phi->vals[5] = state3.Phi->vals[5] = 1;
   state1.Phi->vals[8] = state2.Phi->vals[8] = state3.Phi->vals[8] = 1;
   state1.P->vals[0] = 100;
   state2.P->vals[0] = 0.1;
   state3.P->vals[0] = 0.0;
   state1.P->vals[4] = 100;
   state2.P->vals[4] = 0.1;
   state3.P->vals[4] = 0.0;
   state1.P->vals[8] = 100;
   state2.P->vals[8] = 0.1;
   state3.P->vals[8] = 0.0;

   FILE * f_3 = fopen("dat/second_kp11_p100.dat","w");
   FILE * f_4 = fopen("dat/second_kp11_p0.1.dat","w");
   FILE * f_5 = fopen("dat/second_kp11_p0.dat","w");
   FILE * f_6 = fopen("dat/second_kp22_p100.dat","w");
   FILE * f_7 = fopen("dat/second_kp22_p0.1.dat","w");
   FILE * f_8 = fopen("dat/second_kp22_p0.dat","w");
   FILE * f_9 = fopen("dat/second_kp33_p100.dat","w");
   FILE * f_a = fopen("dat/second_kp33_p0.1.dat","w");
   FILE * f_b = fopen("dat/second_kp33_p0.dat","w");

   for (int i = 1; i < 100; ++i) {

      step_second_order_polynomial_kalman_filter(1.0,1.0,&state1);
      step_second_order_polynomial_kalman_filter(1.0,1.0,&state2);
      step_second_order_polynomial_kalman_filter(1.0,1.0,&state3);
      fprintf(f_3,"%lf  %lf\n",(double)i,state1.P->vals[0]);
      fprintf(f_6,"%lf  %lf\n",(double)i,state1.P->vals[4]);
      fprintf(f_9,"%lf  %lf\n",(double)i,state1.P->vals[8]);
      fprintf(f_4,"%lf  %lf\n",(double)i,state2.P->vals[0]);
      fprintf(f_7,"%lf  %lf\n",(double)i,state2.P->vals[4]);
      fprintf(f_a,"%lf  %lf\n",(double)i,state2.P->vals[8]);
      fprintf(f_5,"%lf  %lf\n",(double)i,state3.P->vals[0]);
      fprintf(f_8,"%lf  %lf\n",(double)i,state3.P->vals[4]);
      fprintf(f_b,"%lf  %lf\n",(double)i,state3.P->vals[8]);
   }

   fclose(f_3);
   fclose(f_4);
   fclose(f_5);
   fclose(f_6);
   fclose(f_7);
   fclose(f_8);
   fclose(f_9);
   fclose(f_a);
   fclose(f_b);

   destroy_second_order_polynomial_kalman_filter(&state1);
   destroy_second_order_polynomial_kalman_filter(&state2);
   destroy_second_order_polynomial_kalman_filter(&state3);

   /* now try and track the falling object xD */
   sopkfs state4;
   init_second_order_polynomial_kalman_filter(&state4);

   FILE * c_d = fopen("dat/second_obj_m.dat","w");
   FILE * c_0 = fopen("dat/second_obj_p.dat","w");
   FILE * c_1 = fopen("dat/second_obj_v.dat","w");
   FILE * c_2 = fopen("dat/second_obj_a.dat","w");

   double ts = 0.1;
   double sigma = 1000.0;
   state4.P->vals[0] = 999999999999999;
   state4.P->vals[4] = 999999999999999;
   state4.P->vals[8] = 999999999999999;
   state4.Phi->vals[0] = 1.0;
   state4.Phi->vals[1] = ts;
   state4.Phi->vals[2] = 0.5*ts*ts;
   state4.Phi->vals[4] = 1.0;
   state4.Phi->vals[5] = ts;
   state4.Phi->vals[8] = 1.0;
   state4.H->vals[0] = 1.0;
   state4.R = sigma*sigma;
   state4.Q = 0.0;

   for (int i = 1; i <= 300; ++i) {

      double t = i*ts;
      double x = 400000 - 6000*t - (32.2*t*t)/2.0;
      double xn = x + grand(0.0,sigma);
      step_second_order_polynomial_kalman_filter(xn,ts,&state4);

      fprintf(c_d,"%lf  %lf\n",t,xn);
      fprintf(c_0,"%lf  %lf\n",t,state4.A->vals[0]);
      fprintf(c_1,"%lf  %lf\n",t,state4.A->vals[1]);
      fprintf(c_2,"%lf  %lf\n",t,state4.A->vals[2]);

   }

   fclose(c_d);
   fclose(c_0);
   fclose(c_1);
   fclose(c_2);

   destroy_second_order_polynomial_kalman_filter(&state4);

   return 0;
}
