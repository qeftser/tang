
/* implimentation of the final fading memory filter presented in 
 * chapter 15. This filter is of interest to me for it's ability
 * to perform relativly well while presenting much less of a 
 * computational burden. This could be put to use on UAVs that
 * have low compute capabilities.
 */

#include "../numerical_basics.h"
#include <stdio.h>
#include <math.h>

DEFMAT(2,1);
DEFMAT(1,2);
DEFMAT(2,2);

int main(void) {

   SRAND();

   double t = 0.0;
   double ts = 0.1;
   double a0 = 400000;
   double a1 = -6000;
   double a2 = -16.1;

   double x = a0 + a1 * t + a2 * t*t;
   double xd = a1 + 2*a2*t;
   double beta = 0.99;
   double gfil = 1.0 - beta;
   double hfil = (1.0 - beta)*(1.0 - beta);
   double xh = 0;
   double xdh = 0;
   double signoise = 1000.0;

   int qswitch = 1;

   double s = 0;

   mat2x1 K = MAT(2,1);
   K.vals[0][0] = gfil;
   K.vals[1][0] = hfil / ts;

   mat1x2 KT = MAT(1,2);
   mattrn((MB *)&K,(MB *)&KT);

   mat2x2 Phi  = MAT(2,2);
   mat2x2 PhiT = MAT(2,2);
   mat2x2 I    = MAT(2,2);
   mat2x2 PhiP = MAT(2,2);
   mat2x2 M    = MAT(2,2);
   mat2x2 KH   = MAT(2,2);
   mat2x2 IKHT = MAT(2,2);
   mat2x2 IKHM = MAT(2,2);
   mat2x2 PRIM = MAT(2,2);
   double R;
   mat2x1 KR   = MAT(2,1);
   mat2x2 KRKT = MAT(2,2);
   mat2x2 P    = MAT(2,2);
   mat1x2 H    = MAT(1,2);

   matzro((MB *)&Phi);
   matzro((MB *)&I);
   matzro((MB *)&M);
   matzro((MB *)&P);
   matzro((MB *)&H);

   I.vals[0][0] = 1.0;
   I.vals[1][1] = 1.0;
   P.vals[0][0] = 9999999999;
   P.vals[1][1] = 9999999999;
   Phi.vals[0][0] = 1.0;
   Phi.vals[0][1] = ts;
   Phi.vals[1][1] = 1.0;
   H.vals[0][0] = 1.0;
   R = signoise*signoise;
   mattrn((MB *)&Phi,(MB *)&PhiT);

   int xn = 0;

   FILE * f_0 = fopen("dat/fading.dat","w");

   while (t < 300.0) {
      xn = xn+1;

      double xk1 = 2.0*(2.0*xn-1)/(xn*(xn+1.0));
      double xk2 = 6.0/(xn*(xn+1.0)*ts);
      matmul((MB *)&Phi,(MB *)&P,(MB *)&PhiP);
      matmul((MB *)&PhiP,(MB *)&PhiT,(MB *)&M);
      matmul((MB *)&K,(MB *)&H,(MB *)&KH);
      matsub((MB *)&I,(MB *)&KH,(MB *)&KH);
      mattrn((MB *)&KH,(MB *)&IKHT);
      matmul((MB *)&KH,(MB *)&M,(MB *)&IKHM);
      matmul((MB *)&IKHM,(MB *)&IKHT,(MB *)&PRIM);
      matsca((MB *)&K,(MB *)&KR,R);
      matmul((MB *)&KR,(MB *)&KT,(MB *)&KRKT);
      matadd((MB *)&PRIM,(MB *)&KRKT,(MB *)&P);

      double xnoise = grand(0.0,signoise);

      x = a0+a1*t+a2*t*t;
      xd = a1+2*a2*t;
      double xs = x+xnoise;

      double res = xs-xh-ts*xdh;

      double gain1, gain2, sp00, sp11;
      if (xk1 > gfil) {
         gain1 = xk1;
         gain2 = xk2;
         if (xn == 1.0) {
            sp00 = 0.0;
            sp11 = 0.0;
         }
         else {
            sp00 = signoise*sqrt(2.0*(2.0*xn-1)/(xn*(xn+1.0)));
            sp11 = signoise*sqrt(12/(xn*(xn*xn-1)*ts*ts));
         }
      }
      else {
         if (qswitch) {
            qswitch = 0.0;
            P.vals[0][0] = sp00*sp00;
            P.vals[1][1] = sp11*sp11;
            P.vals[0][1] = P.vals[1][0] = 0.0;
         }
         gain1 = gfil;
         gain2 = hfil/ts;
         sp00 = sqrt(P.vals[0][0]);
         sp11 = sqrt(P.vals[1][1]);
      }

      xh = xh+xdh*ts-16.1*ts*ts+gain1*res;
      xdh = xdh-32.2*ts+gain2*res;

      double xherr = x-xh;
      double xdherr =xd-xdh;

      fprintf(f_0,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
              t,x,xd,xs,xh,xdh,xherr,xdherr,sp00,-sp00,sp11,-sp11);

      t += ts;
   }

   fclose(f_0);

   return 0;
}
