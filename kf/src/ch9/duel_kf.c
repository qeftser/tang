
/* implimentation of the decoupled kalman filters
 * presented in chapter 9. These have a cool
 * name so I hope they are cool to impliment
 */

#include "../numerical_basics.h"
#include <stdio.h>
#include <math.h>

DEFMAT(2,2);
DEFMAT(1,2);
DEFMAT(2,1);
DEFMAT(1,1);

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
   double xr = 100000.0;
   double yr = 0.0;
   double t = 0.0;
   double s = 0.0;
   double h = 0.001;

   mat2x2 Phi = MAT(2,2);
   matzro((MB *)&Phi);
   mat2x2 PX = MAT(2,2);
   matzro((MB *)&PX);
   mat2x2 PY = MAT(2,2);
   matzro((MB *)&PY);
   mat2x2 Q = MAT(2,2);
   matzro((MB *)&Q);
   mat2x2 I = MAT(2,2);
   matzro((MB *)&I);
   mat1x2 H = MAT(1,2);

   mat2x2 MX = MAT(2,2);
   mat2x2 MY = MAT(2,2);
   mat2x1 KX = MAT(2,1);
   mat2x1 KY = MAT(2,1);

   Phi.vals[0][0] = 1.0;
   Phi.vals[0][1] = ts;
   Phi.vals[1][1] = 1.0;
   H.vals[0][0] = 1.0;
   H.vals[0][1] = 0.0;

   mat2x2 PhiT = MAT(2,2);
   mat2x1 HT = MAT(2,1);
   mattrn((MB *)&Phi,(MB *)&PhiT);
   mattrn((MB *)&H,(MB *)&HT);

   I.vals[0][0] = 1.0;
   I.vals[1][1] = 1.0;
   Q.vals[0][0] = pnoise*ts*ts*ts/3.0;
   Q.vals[0][1] = pnoise*ts*ts/2.0;
   Q.vals[1][0] = Q.vals[0][1];
   Q.vals[1][1] = pnoise*ts;
   PX.vals[0][0] = 1000.0 * 1000.0;
   PX.vals[1][1] = 100.0 * 100.0;
   PY.vals[0][0] = 1000.0 * 1000.0;
   PY.vals[1][1] = 100.0 * 100.0;
   
   double xth = xt + 1000;
   double xtdh = xtd - 100;
   double yth = yt - 1000;
   double ytdh = ytd + 100;

   FILE * f_0 = fopen("dat/duel_kf.dat","w");

   double xtold, xtdold, ytold, ytdold, xtdd, ytdd;
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
      xt = 0.5 * (xtold+xt+h*xtd);
      xtd = 0.5 * (xtdold+xtd+h*xtdd);
      yt = 0.5 * (ytold+yt+h*ytd);
      ytd = 0.5 * (ytdold+ytd+h*ytdd);
      s = s+h;

      if (s >= (ts-0.00001)) {
         s = 0.0;

         double theth = atan2(yth-yr,xth-xr);
         double rth = sqrt((xth-xr)*(xth-xr)+(yth-yr)*(yth-yr));
         double R = (cos(theth)*sigr)*(cos(theth)*sigr)+(rth*sin(theth)*sigth)*(rth*sin(theth)*sigth);

         mat2x2 PhiP = MAT(2,2);
         mat2x2 PhiPPhiT = MAT(2,2);
         matmul((MB *)&Phi,(MB *)&PX,(MB *)&PhiP);
         matmul((MB *)&PhiP,(MB *)&PhiT,(MB *)&PhiPPhiT);
         matadd((MB *)&PhiPPhiT,(MB *)&Q,(MB *)&MX);

         mat1x2 HM = MAT(1,2);
         mat1x1 HMHT = MAT(1,1);
         mat2x1 MHT = MAT(2,1);
         matmul((MB *)&H,(MB *)&MX,(MB *)&HM);
         matmul((MB *)&HM,(MB *)&HT,(MB *)&HMHT);
         HMHT.vals[0][0] += R;
         matmul((MB *)&MX,(MB *)&HT,(MB *)&MHT);
         matsca((MB *)&MHT,(MB *)&KX,1.0 / HMHT.vals[0][0]);

         mat2x2 KH = MAT(2,2);
         matmul((MB *)&KX,(MB *)&H,(MB *)&KH);
         matsub((MB *)&I,(MB *)&KH,(MB *)&KH);
         matmul((MB *)&KH,(MB *)&MX,(MB *)&PX);

         double thet = atan2(yt-yr,xt-xr);
         double rt = sqrt((xt-xr)*(xt-xr)+(yt-yr)*(yt-yr));
         double thetmeas = thet + grand(0.0,sigth);
         double rtmeas = rt + grand(0.0,sigr);
         double xtmeas = rtmeas*cos(thetmeas)+xr;
         double res1 = xtmeas-xth-ts*xtdh;

         xth = xth+ts*xtdh+KX.vals[0][0]*res1;
         xtdh = xtdh+KX.vals[1][0]*res1;

         R = (sin(theth)*sigr)*(sin(theth)*sigr)+(rth*cos(theth)*sigth)*(rth*cos(theth)*sigth);
         
         matmul((MB *)&Phi,(MB *)&PY,(MB *)&PhiP);
         matmul((MB *)&PhiP,(MB *)&PhiT,(MB *)&PhiPPhiT);
         matadd((MB *)&PhiPPhiT,(MB *)&Q,(MB *)&MY);

         matmul((MB *)&H,(MB *)&MY,(MB *)&HM);
         matmul((MB *)&HM,(MB *)&HT,(MB *)&HMHT);
         HMHT.vals[0][0] += R;
         matmul((MB *)&MY,(MB *)&HT,(MB *)&MHT);
         matsca((MB *)&MHT,(MB *)&KY,1.0 / HMHT.vals[0][0]);

         matmul((MB *)&KY,(MB *)&H,(MB *)&KH);
         matsub((MB *)&I,(MB *)&KH,(MB *)&KH);
         matmul((MB *)&KH,(MB *)&MY,(MB *)&PY);

         double ytmeas = rtmeas*sin(thetmeas)+yr;
         double res2 = ytmeas-yth-ts*ytdh+0.5*ts*ts*g;
         
         yth = yth+ts*ytdh-0.5*ts*ts*g+KY.vals[0][0]*res2;
         ytdh = ytdh-ts*g+KY.vals[1][0]*res2;

         double SP00 = sqrt(PX.vals[0][0]);
         double SP11 = sqrt(PX.vals[1][1]);
         double SP22 = sqrt(PY.vals[0][0]);
         double SP33 = sqrt(PY.vals[1][1]);

         fprintf(f_0,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                 t, xt, xtd, yt, ytd, thet*57.3, rt, xth, xtdh, yth, ytdh, xth-xt, xtdh-xtd, yth-yt, ytdh-ytd,
                 SP00,-SP00,SP11,-SP11,SP22,-SP22,SP33,-SP33);

      }

   }

   fclose(f_0);

   return 0;
}
