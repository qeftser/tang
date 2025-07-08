
/* implimentation of the extended kalman filter 
 * presented in chapter 9. This is the one in
 * cartesian coordinates. It is going to be copied
 * from the source, as I am tired of introducing 
 * small errors that I can't fix.
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

   double ts = 1.0;
   double pnoise = 0.0;
   double sigth = 0.01;
   double sigr = 100.0;
   double vt = 3000.0;
   double gamdeg = 45.0;
   double g = 32.2;
   double xt = 0.0;
   double yt = 0.0;
   double xtd = vt*cos(gamdeg/57.3);
   double ytd = vt*sin(gamdeg/57.3);
   double xr = 100000;
   double yr = 0.0;
   double t = 0.0;
   double s = 0.0;
   double h = 0.001;

   mat4x4 Phi = MAT(4,4);
   mat4x4 P   = MAT(4,4);
   mat4x4 Q   = MAT(4,4);
   mat4x4 I   = MAT(4,4);
   mat2x4 H   = MAT(2,4);
   mat4x2 HT  = MAT(4,2);
   mat4x4 M   = MAT(4,4);
   mat4x2 K   = MAT(4,2);

   mat4x4 PhixP = MAT(4,4);
   mat4x4 PhixPxPhiT = MAT(4,4);
   mat2x4 HxM = MAT(2,4);
   mat2x2 HxMxHT = MAT(2,2);
   mat2x2 HxMxHTpR_INV = MAT(2,2);
   mat4x2 MxHT = MAT(4,2);
   mat4x4 KxH  = MAT(4,4);

   double ts2 = ts*ts;
   double ts3 = ts*ts*ts;

   Q.vals[0][0] = pnoise*ts3/3.0;
   Q.vals[0][1] = pnoise*ts2/2.0;
   Q.vals[1][0] = Q.vals[0][1];
   Q.vals[1][1] = pnoise*ts;
   Q.vals[2][2] = Q.vals[0][0];
   Q.vals[2][3] = Q.vals[0][1];
   Q.vals[3][2] = Q.vals[1][0];
   Q.vals[3][3] = Q.vals[1][1];

   matzro((MB *)&Phi);
   Phi.vals[0][0] = 1.0;
   Phi.vals[0][1] = ts;
   Phi.vals[1][1] = 1.0;
   Phi.vals[2][2] = 1.0;
   Phi.vals[2][3] = ts;
   Phi.vals[3][3] = 1.0;

   mat4x4 PhiT = MAT(4,4);
   mattrn((MB *)&Phi,(MB *)&PhiT);

   mat2x2 R = MAT(2,2);
   R.vals[0][0] = sigth*sigth;
   R.vals[1][0] = R.vals[0][1] = 0.0;
   R.vals[1][1] = sigr*sigr;

   matzro((MB *)&I);
   I.vals[0][0] = 1.0;
   I.vals[1][1] = 1.0;
   I.vals[2][2] = 1.0;
   I.vals[3][3] = 1.0;

   matzro((MB *)&P);
   P.vals[0][0] = 1000.0*1000.0;
   P.vals[1][1] = 100.0*100.0;
   P.vals[2][2] = 1000.0*1000.0;
   P.vals[3][3] = 100.0*100.0;

   double xth = xt + 1000;
   double xtdh = xtd - 100;
   double yth = yt - 1000;
   double ytdh = ytd + 100;

   FILE * f_0 = fopen("dat/ekf.dat","w");

   double xtold, xtdold, ytold, ytdold, xtdd, ytdd, xtb, ytb, xtdb, ytdb, rtb;
   while (yt >= 0.0) {

      xtold = xt;
      xtdold = xtd;
      ytold = yt;
      ytdold = ytd;
      xtdd = 0.0;
      ytdd = -g;

      xt = xt+h*xtd;
      xtd = xtd+h*xtdd;
      yt = yt+h*ytd;
      ytd = ytd+h*ytdd;

      t = t+h;

      xt = 0.5 * (xtold + xt+h*xtd);
      xtd = 0.5 * (xtdold + xtd+h*xtdd);
      yt = 0.5 * (ytold + yt+h*ytd);
      ytd = 0.5 * (ytdold + ytd+h*ytdd);

      s = s+h;

      if (s >= (ts - 0.00001)) {
         s = 0.0;

         xtb = xth + ts*xtdh;
         xtdb = xtdh;
         ytb = yth + ts*ytdh - 0.5*g*ts*ts;
         ytdb = ytdh - g*ts;
         rtb = sqrt((xtb-xr)*(xtb-xr) + (ytb-yr)*(ytb-yr));

         H.vals[0][0] = -(ytb-yr)/(rtb*rtb);
         H.vals[0][1] = 0.0;
         H.vals[0][2] = (xtb-xr)/(rtb*rtb);
         H.vals[0][3] = 0.0;
         H.vals[1][0] = (xtb-xr)/rtb;
         H.vals[1][1] = 0.0;
         H.vals[1][2] = (ytb-yr)/rtb;
         H.vals[1][3] = 0.0;
         mattrn((MB *)&H,(MB *)&HT);

         matmul((MB *)&Phi,(MB *)&P,(MB *)&PhixP);
         matmul((MB *)&PhixP,(MB *)&PhiT,(MB *)&PhixPxPhiT);
         matadd((MB *)&PhixPxPhiT,(MB *)&Q,(MB *)&M);

         matmul((MB *)&H,(MB *)&M,(MB *)&HxM);
         matmul((MB *)&HxM,(MB *)&HT,(MB *)&HxMxHT);
         matadd((MB *)&HxMxHT,(MB *)&R,(MB *)&HxMxHT);
         matinv((MB *)&HxMxHT,(MB *)&HxMxHTpR_INV);
         matmul((MB *)&M,(MB *)&HT,(MB *)&MxHT);
         matmul((MB *)&MxHT,(MB *)&HxMxHTpR_INV,(MB *)&K);

         matmul((MB *)&K,(MB *)&H,(MB *)&KxH);
         matsub((MB *)&I,(MB *)&KxH,(MB *)&KxH);
         matmul((MB *)&KxH,(MB *)&M,(MB *)&P);

         double th = atan2(yt-yr,xt-xr);
         double r = sqrt((xt-xr)*(xt-xr) + (yt-yr)*(yt-yr));

         double zth = th + grand(0.0,sigth);
         double zr  = r  + grand(0.0,sigr);

         double thb = atan2(ytb-yr,xtb-xr);
         double rb  = sqrt((xtb-xr)*(xtb-xr) + (ytb-yr)*(ytb-yr));

         double res1 = zth - thb;
         double res2 = zr  - rb;

         xth  = xtb  + K.vals[0][0]*res1 + K.vals[0][1]*res2;
         xtdh = xtdb + K.vals[1][0]*res1 + K.vals[1][1]*res2;
         yth  = ytb  + K.vals[2][0]*res1 + K.vals[2][1]*res2;
         ytdh = ytdb + K.vals[3][0]*res1 + K.vals[3][1]*res2;

         double SP00 = sqrt(P.vals[0][0]);
         double SP11 = sqrt(P.vals[1][1]);
         double SP22 = sqrt(P.vals[2][2]);
         double SP33 = sqrt(P.vals[3][3]);

         fprintf(f_0,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                 t, xt, xtd, yt, ytd, th*57.3, r, xth, xtdh, yth, ytdh, xth-xt, xtdh-xtd, yth-yt, ytdh-ytd,
                 SP00,-SP00,SP11,-SP11,SP22,-SP22,SP33,-SP33);

      }
   }

   fclose(f_0);

   return 0;
}
