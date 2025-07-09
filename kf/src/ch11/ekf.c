
/* reproduction of the extended kalman filter estimating the
 * position of an object moving with a variable velocity
 * using the measurements of two satellites. It is pretty
 * easy to extend the filter and add or remove measurements
 * once a basic version is completed. Like all things, I guess
 * it is easier to start small and expand then do it all at once...
 */

#include "../numerical_basics.h"
#include <stdio.h>
#include <math.h>

DEFMAT(4,4);
DEFMAT(4,2);
DEFMAT(2,4);
DEFMAT(2,2);

int main(void) {

   SRAND();

   double signoise = 300.0;
   double phireal = 9000.0;
   double tau = 5.0;
   double phis = 100.0;
   double x = 0.0;
   double y = 0.0;
   double y1 = 0.0;
   double xdp = 100.0;
   double yd = 0.0;
   double xh = 1000.0;
   double yh = 2000.0;
   double xdh = 0.0;
   double ydh = 0.0;
   double xr1 = 1000000;
   double yr1 = 20000 * 3280.0;
   double xr2 = 500000000;
   double yr2 = 20000 * 3280.0;
   double xr2d = -14600;
   double xr1d = -14600;
   double rs = 1.0;
   double tf = 200;
   double t = 0.0;
   double s = 0.0;
   double h = 0.01;
   double sig = sqrt(phireal/h);
   double ts = 1.0;
   double ts2 = ts*ts;
   double ts3 = ts*ts*ts;

   mat4x4 Phi  = MAT(4,4);
   mat4x4 P    = MAT(4,4);
   mat4x4 M    = MAT(4,4);
   mat4x4 Q    = MAT(4,4);
   mat4x4 I    = MAT(4,4);
   mat4x4 PhiT = MAT(4,4);
   mat2x4 H    = MAT(2,4);
   mat4x2 HT   = MAT(4,2);
   mat4x2 K    = MAT(4,2);
   mat2x2 R    = MAT(2,2);

   matzro((MB *)&Phi);
   matzro((MB *)&P);
   matzro((MB *)&M);
   matzro((MB *)&Q);
   matzro((MB *)&H);
   matzro((MB *)&K);
   matzro((MB *)&I);
   matzro((MB *)&R);

   I.vals[0][0] = 1.0;
   I.vals[1][1] = 1.0;
   I.vals[2][2] = 1.0;
   I.vals[3][3] = 1.0;
   P.vals[0][0] = 1000.0 * 1000.0;
   P.vals[1][1] = 100.0 * 100.0;
   P.vals[2][2] = 2000.0 * 2000.0;
   P.vals[3][3] = 100.0 * 100.0;
   Q.vals[0][0] = phis*ts3/3.0;
   Q.vals[0][1] = phis*ts2/2.0;
   Q.vals[1][0] = Q.vals[0][1];
   Q.vals[1][1] = phis*ts;
   Q.vals[2][2] = Q.vals[0][0];
   Q.vals[2][3] = Q.vals[0][1];
   Q.vals[3][2] = Q.vals[1][0];
   Q.vals[3][3] = Q.vals[1][1];
   R.vals[0][0] = signoise*signoise;
   R.vals[0][1] = 0.0;
   R.vals[1][0] = 0.0;
   R.vals[1][1] = signoise*signoise;

   FILE * f_0 = fopen("dat/ekf.dat","w");

   while (t <= tf) {

      double x1 = grand(0.0,sig);
      double xr1old = xr1;
      double xr2old = xr2;
      double xold = x;
      double yold = y;
      double y1old = y1;
      xr1d = -14600;
      xr2d = -14600;
      double y1d = (x1-y1)/tau;
      double xd = xdp + y1;

      yd = 0;
      xr1 = xr1+h*xr1d;
      xr2 = xr2+h*xr2d;
      x = x+h*xd;
      y = y+h*yd;
      y1 = y1+h*y1d;
      t = t+h;
      xr1d = -14600;
      xr2d = -14600;
      y1d = (x1-y1)/tau;
      xd = xdp+y1;
      yd = 0;
      xr1 = 0.5*(xr1old+xr1+h*xr1d);
      xr2 = 0.5*(xr2old+xr2+h*xr2d);
      x = 0.5*(xold+x+h*xd);
      y = 0.5*(yold+y+h*yd);
      y1 = 0.5*(y1old+y1+h+y1d);
      s = s+h;

      if (s >= (ts-0.00001)) {
         s = 0.0;

         double xb = xh+xdh*ts;
         double yb = yh+ydh*ts;
         double r1b = sqrt((xr1-xb)*(xr1-xb)+(yr1-yb)*(yr1-yb));
         double r2b = sqrt((xr2-xb)*(xr2-xb)+(yr2-yb)*(yr2-yb));
         H.vals[0][0] = -(xr1-xb)/r1b;
         H.vals[0][1] = 0.0;
         H.vals[0][2] = -(yr1-yb)/r1b;
         H.vals[0][3] = 0.0;
         H.vals[1][0] = -(xr2-xb)/r2b;
         H.vals[1][1] = 0.0;
         H.vals[1][2] = -(yr2-yb)/r2b;
         H.vals[1][3] = 0.0;
         mattrn((MB *)&H,(MB *)&HT);

         Phi.vals[0][0] = 1.0;
         Phi.vals[0][1] = ts;
         Phi.vals[1][1] = 1.0;
         Phi.vals[2][2] = 1.0;
         Phi.vals[2][3] = ts;
         Phi.vals[3][3] = 1.0;
         mattrn((MB *)&Phi,(MB *)&PhiT);

         mat4x4 PhiP = MAT(4,4);
         mat4x4 PhiPPhiT = MAT(4,4);
         mat2x4 HM = MAT(2,4);
         mat2x2 HMHT = MAT(2,2);
         mat2x2 HMHTRINV = MAT(2,2);
         mat4x2 MHT = MAT(4,2);
         mat4x4 KH = MAT(4,4);

         matmul((MB *)&Phi,(MB *)&P,(MB *)&PhiP);
         matmul((MB *)&PhiP,(MB *)&PhiT,(MB *)&PhiPPhiT);
         matadd((MB *)&PhiPPhiT,(MB *)&Q,(MB *)&M);

         matmul((MB *)&H,(MB *)&M,(MB *)&HM);
         matmul((MB *)&HM,(MB *)&HT,(MB *)&HMHT);
         matadd((MB *)&HMHT,(MB *)&R,(MB *)&HMHT);
         matinv((MB *)&HMHT,(MB *)&HMHTRINV);
         matmul((MB *)&M,(MB *)&HT,(MB *)&MHT);
         matmul((MB *)&MHT,(MB *)&HMHTRINV,(MB *)&K);

         matmul((MB *)&K,(MB *)&H,(MB *)&KH);
         matsub((MB *)&I,(MB *)&KH,(MB *)&KH);
         matmul((MB *)&KH,(MB *)&M,(MB *)&P);

         double r1noise = grand(0.0,signoise);
         double r2noise = grand(0.0,signoise);
         double r1 = sqrt((xr1-x)*(xr1-x)+(yr1-y)*(yr1-y));
         double r2 = sqrt((xr2-x)*(xr2-x)+(yr2-y)*(yr2-y));
         double res1 = r1 + r1noise-r1b;
         double res2 = r2 + r2noise-r2b;

         xh = xb + K.vals[0][0]*res1 + K.vals[0][1]*res2;
         xdh = xdh + K.vals[1][0]*res1 + K.vals[1][1]*res2;
         yh = yb + K.vals[2][0]*res1 + K.vals[2][1]*res2;
         ydh = ydh + K.vals[3][0]*res1 + K.vals[3][1]*res2;

         double errx = x-xh;
         double sp00 = sqrt(P.vals[0][0]);
         double errxd = xd-xdh;
         double sp11 = sqrt(P.vals[1][1]);
         double erry = y-yh;
         double sp22 = sqrt(P.vals[2][2]);
         double erryd = yd-ydh;
         double sp33 = sqrt(P.vals[3][3]);

         fprintf(f_0,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                 t,x,y,xd,yd,xh,yh,xdh,ydh,errx,erry,errxd,erryd,sp00,-sp00,sp11,-sp11,sp22,-sp22,sp33,-sp33);

      }

   }

   fclose(f_0);

   return 0;
}

