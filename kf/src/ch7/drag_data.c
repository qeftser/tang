
/* implimentation of the extended kalman filters described in chapter
 * 7. Here things are kind of different in that the ekf functions
 * f and h are custom and need to be recomputed at every stage. Because
 * of this, we cannot really encapsulate things into a general class,
 * so we will do everything in this file here :). I promise I am not
 * being lazy! Please!
 */

#include "../numerical_basics.h"
#include <stdio.h>
#include <math.h>

DEFMAT(2,2);
DEFMAT(2,1);
DEFMAT(1,2);
DEFMAT(1,1);

void project(double ts, double xp, double xdp, double beta, double * xh, double * xdh, double * xddh, double hp) {

   double t = 0;
   double x = xp;
   double xd = xdp;
   double h = hp;
   double xdd;

   while (t <= (ts - 0.0001)) {
      xdd = 0.0034*32.2*xd*xd*exp(-x/22000.0)/(2.0*beta)-32.2;
      xd = xd + h*xdd;
      x = x+h*xd;

      t += h;
   }

   *xh = x;
   *xdh = xd;
   *xddh = xdd;
}

int main(void) {

   SRAND();

   double _infinity = 999999999999999;
   double sigma = 1000;
   double x, xd, xdd, beta, ts, tf, t, s, h, x_old, xd_old, g;

   double pnoise = 0.0;

   /* NOTE: here changing the beta term changes the drag
    *       on the falling object. */
   beta = 500;

   ts = 0.1;
   tf = 40;
   h = 0.001;
   g = 32.2;

   mat2x2 P = MAT(2,2);
   mat2x1 K = MAT(2,1);
   mat2x2 M = MAT(2,2);

   P.vals[0][0] = _infinity;
   P.vals[1][1] = _infinity;

   double a0 = 0.0, a1 = 0.0;

   FILE * f_0 = fopen("dat/object_falling.dat","w");
   FILE * f_1 = fopen("dat/object_tracking.dat","w");

   t = 0; s = 0; x = 200000, xd = -6000;
   while (t < tf) {

      x_old = x;
      xd_old = xd;
      xdd = 0.0034*g*xd*xd*exp(-x/22000.0)/(2.0*beta)-g;
      x = x + h*xd;
      xd = xd + h*xdd;

      t += h;

      xdd = 0.0034*g*xd*xd*exp(-x/22000.0)/(2.0*beta)-g;
      x = 0.5 * (x_old+x+h*xd);
      xd = 0.5 * (xd_old+xd+h*xdd);

      s += h;

      if (s >= ts - 0.00001) {
         s = 0;

         double z = x + grand(0.0,sigma);

         fprintf(f_0,"%lf  %lf  %lf  %lf\n",t,x,xd,xdd);

         double f21 = (-0.0034*g*a1*a1*exp(-a0/22000.0)/(44000*beta));
         double f22 = (a1*g*0.0034*exp(-a0/22000.0)/beta);
         double R = sigma*sigma;

         /* compute the matrices for this timestep */

         mat2x2 Phi = MAT(2,2);
         Phi.vals[0][0] = 1.0;
         Phi.vals[0][1] = ts;
         Phi.vals[1][0] = f21 * ts;
         Phi.vals[1][1] = f22 * ts + 1.0;

         mat1x2 H = MAT(1,2);
         H.vals[0][0] = 1.0;
         H.vals[0][1] = 0.0;

         mat2x2 Q = MAT(2,2);
         Q.vals[0][0] = ts*ts*ts / 3.0;
         Q.vals[0][1] = (ts*ts / 2.0) + f22*(ts*ts*ts / 3.0);
         Q.vals[1][0] = (ts*ts / 2.0) + f22*(ts*ts*ts / 3.0);
         Q.vals[1][1] = ts + f22*(ts*ts) + f22*f22*(ts*ts*ts / 3.0);
         matsca((matbase *)&Q,(matbase *)&Q,pnoise);

         mat2x2 PhiT = MAT(2,2);
         mat2x1 HT = MAT(2,1);
         mattrn((matbase *)&Phi,(matbase *)&PhiT);
         mattrn((matbase *)&H,(matbase *)&HT);

         mat2x2 interm0 = MAT(2,2);
         mat2x1 interm1 = MAT(2,1);
         mat1x1 interm2 = MAT(1,1);

         /* compute the matrix riccati equations */
         matmul((matbase *)&Phi,(matbase *)&P,(matbase *)&interm0);
         matmul((matbase *)&interm0,(matbase *)&PhiT,(matbase *)&M);
         matadd((matbase *)&Q,(matbase *)&M,(matbase *)&M);

         matmul((matbase *)&M,(matbase *)&HT,(matbase *)&K);
         matmul((matbase *)&M,(matbase *)&HT,(matbase *)&interm1);
         matmul((matbase *)&H,(matbase *)&interm1,(matbase *)&interm2);
         matsca((matbase *)&K,(matbase *)&K, 1.0 / (interm2.vals[0][0] + R));

         matmul((matbase *)&K,(matbase *)&H,(matbase *)&interm0);
         interm0.vals[0][0] = 1.0 - interm0.vals[0][0];
         interm0.vals[0][1] = - interm0.vals[0][1];
         interm0.vals[1][0] = - interm0.vals[1][0];
         interm0.vals[1][1] = 1.0 - interm0.vals[1][1];
         matmul((matbase *)&interm0,(matbase *)&M,(matbase *)&P);

         /* update the states given this information */
         double xdd_bar = 0.0034*g*a1*a1*exp(-a0/22000.0)/(2.0*beta) - g;
         double xd_bar;
         double x_bar;

         project(ts,a0,a1,beta,&x_bar,&xd_bar,&xdd_bar,h);

         a0 = x_bar + K.vals[0][0] * (z - x_bar);
         a1 = xd_bar + K.vals[1][0] * (z - x_bar);

         fprintf(f_1,"%lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n",t,a0,a0-x,a1,a1-xd,
                 sqrt(P.vals[0][0]),-sqrt(P.vals[0][0]),sqrt(P.vals[1][1]),-sqrt(P.vals[1][1]));

      }

   }

   fclose(f_0);
   fclose(f_1);

   return 0;
}
