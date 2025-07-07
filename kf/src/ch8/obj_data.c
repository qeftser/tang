
/* Implimentation of the functional nonlinear kalman filter 
 * presented in chapter 8. */

#include "../numerical_basics.h"
#include <stdio.h>
#include <math.h>

DEFMAT(3,3);
DEFMAT(3,1);
DEFMAT(1,3);
DEFMAT(1,1);

void project(double x, double xd, double beta, double ts, double h, double * x_bar, double * xd_bar, double * xdd_bar) {

   double t = 0.0;
   double xdd, x_old, xd_old;

   while (t < ts-h/10) {

      x_old = x;
      xd_old = xd;

      xdd = (0.0034*32.2*exp(-x/22000.0)*xd*xd) / (2.0 * beta) - 32.2;
      x = x + xd*h;
      xd = xd + xdd*h;

      t += h;

      xdd = (0.0034*32.2*exp(-x/22000.0)*xd*xd) / (2.0 * beta) - 32.2;
      x = 0.5 * (x_old + x + h*xdd);
      xd = 0.5 * (xd_old + xd + h*xdd);
   }

   *x_bar = x;
   *xd_bar = xd;
   *xdd_bar = xdd;
}

int main(void) {

   SRAND();

   double sigma = 25;
   double _infinity = 999999999999999;
   double pnoise = 1.0;
   double x, xd, xdd, x_bar, xd_bar, xdd_bar, t, ts, h, s, a0, a1, tl, beta, a2, x_old, xd_old;
   double p, f21, f22, f23;

   tl = 30.0;
   t = 0;
   s = 0;
   ts = 0.1;
   h = 0.001;

   a0 = 200025;
   a1 = -6150;
   a2 = 1.0 / 800.0;

   mat1x3 H =    MAT(1,3);
   mat3x3 M =    MAT(3,3);
   mat3x1 K =    MAT(3,1);
   mat3x3 P =    MAT(3,3);
   mat3x3 Phi =  MAT(3,3);
   mat3x3 G =    MAT(3,3);
   mat3x3 Q =    MAT(3,3);
   mat3x3 i0 =   MAT(3,3);
   mat3x1 i1 =   MAT(3,1);
   mat1x1 i2 =   MAT(1,1);
   mat3x1 HT =   MAT(3,1);
   mat3x3 PhiT = MAT(3,3);
   mat3x3 I =    MAT(3,3);
   double R = sigma*sigma;

   matzro((MB *)&P);
   P.vals[0][0] = sigma*sigma;
   P.vals[1][1] = 20000;
   P.vals[2][2] = ((1.0 / 500) - a2) * ((1.0 / 500) - a2);

   matzro((MB *)&I);
   I.vals[0][0] = 1.0;
   I.vals[1][1] = 1.0;
   I.vals[2][2] = 1.0;

   H.vals[0][0] = 1.0;
   H.vals[0][1] = 0.0;
   H.vals[0][2] = 0.0;
   mattrn((MB *)&H,(MB *)&HT);

   x = 200000;
   xd = -6000;
   beta = 500;

   FILE * f_0 = fopen("dat/obj.dat","w");

   while (t < tl) {

      x_old = x;
      xd_old = xd;

      xdd = (0.0034*32.2*exp(-x/22000.0)*xd*xd) / (2.0 * beta) - 32.2;
      x = x + xd*h;
      xd = xd + xdd*h;

      t += h;

      xdd = (0.0034*32.2*exp(-x/22000.0)*xd*xd) / (2.0 * beta) - 32.2;
      x = 0.5 * (x_old + x + h*xdd);
      xd = 0.5 * (xd_old + xd + h*xdd);

      s += h;

      if (s >= (ts - h/10)) {
         s = 0.0;

         double z = x + grand(0.0,sigma);

         project(a0,a1,1.0 / a2,ts,h,&x_bar,&xd_bar,&xdd_bar);

         p = 0.0034*exp(-a0/22000.0);
         f21 = (-p*32.2*a1*a1*a2) / (44000.0);
         f22 = (p*a1*32.2*a2);
         f23 = (p*16.1*a1*a1);

         Phi.vals[0][0] = 1.0;
         Phi.vals[0][1] = ts;
         Phi.vals[0][2] = 0.0;
         Phi.vals[1][0] = f21*ts;
         Phi.vals[1][1] = f22*ts + 1.0;
         Phi.vals[1][2] = f23*ts;
         Phi.vals[2][0] = 0.0;
         Phi.vals[2][1] = 0.0;
         Phi.vals[2][2] = 1.0;
         mattrn((MB *)&Phi,(MB *)&PhiT);

         Q.vals[0][0] = 0.0;
         Q.vals[0][1] = 0.0;
         Q.vals[0][2] = 0.0;
         Q.vals[1][0] = 0.0;
         Q.vals[1][1] = f23*f23*f23*ts*ts*ts/3.0;
         Q.vals[1][2] = f23*ts*ts/2.0;
         Q.vals[2][0] = 0.0;
         Q.vals[2][1] = f23*ts*ts/2.0;
         Q.vals[2][2] = ts;
         matsca((MB *)&Q,(MB *)&Q,pnoise);

         matmul((MB *)&Phi,(MB *)&P,(MB *)&i0);
         matmul((MB *)&i0,(MB *)&PhiT,(MB *)&M);
         matadd((MB *)&M,(MB *)&Q,(MB *)&M);

         matmul((MB *)&M,(MB *)&HT,(MB *)&i1);
         matmul((MB *)&H,(MB *)&i1,(MB *)&i2);
         matsca((MB *)&i1,(MB *)&K,1.0/(i2.vals[0][0]+R));

         matmul((MB *)&K,(MB *)&H,(MB *)&i0);
         matsub((MB *)&I,(MB *)&i0,(MB *)&i0);
         matmul((MB *)&i0,(MB *)&M,(MB *)&P);

         a0 = x_bar + K.vals[0][0] * (z - x_bar);
         a1 = xd_bar + K.vals[1][0] * (z - x_bar);
         a2 = a2 + K.vals[2][0] * (z - x_bar);

         fprintf(f_0,"%lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n",
                 t,x,xd,a0,a1,a2,a0-x,a1-xd,a2-(1.0/beta),sqrt(P.vals[0][0]),-sqrt(P.vals[0][0]),
                 sqrt(P.vals[1][1]),-sqrt(P.vals[1][1]),sqrt(P.vals[2][2]),-sqrt(P.vals[2][2]));

      }
   }

   fclose(f_0);

   return 0;
}
