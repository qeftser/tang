
#include "polynomial_kalman_filters.h"

void step_zeroth_order_polynomial_kalman_filter(double measurement, double elapsed, zopkfs * state) {

   /* first perform the Riccati equations to update the kalman gain */
   {
      state->M = state->Phi * state->P * state->Phi + (state->Q * elapsed);
      state->K = (state->M * state->H) / (state->H * state->M * state->H + state->R);
      state->P = (1 - state->K * state->H) * state->M;
   }

   /* now update our estimate of the state */
   state->a0 = state->Phi * state->a0 + state->K * (measurement - state->H * state->Phi * state->a0);

}

void init_first_order_polynomial_kalman_filter(fopkfs * state) {
   matzro(state->Phi = DALLOCMAT(2,2));
   matzro(state->H = DALLOCMAT(1,2));
   matzro(state->P = DALLOCMAT(2,2));
   matzro(state->M = DALLOCMAT(2,2));
   matzro(state->K = DALLOCMAT(2,1));
   matzro(state->A = DALLOCMAT(2,1));
}

void destroy_first_order_polynomial_kalman_filter(fopkfs * state) {
   free(state->Phi);
   free(state->H);
   free(state->P);
   free(state->M);
   free(state->K);
   free(state->A);
}

void step_first_order_polynomial_kalman_filter(double measurement, double elapsed, fopkfs * state) {

   /* prepare the matrix Qk for this update */
   mat Qk = DALLOCMAT(2,2);
   Qk->vals[0] = (elapsed*elapsed*elapsed) / 3.0;
   Qk->vals[1] = Qk->vals[2] = (elapsed*elapsed) / 2.0;
   Qk->vals[3] = elapsed;
   matsca(Qk,Qk,state->Q);

   {
      mat interm0 = DALLOCMAT(2,2),
          interm1 = DALLOCMAT(2,2),
          interm2 = DALLOCMAT(2,1),
          interm3 = DALLOCMAT(1,2),
          interm4 = DALLOCMAT(1,1),
          interm5 = DALLOCMAT(2,1);

      /* perform the Riccati equations to update the kalman gain */

      /* Mk = Phi * Pk-1 * (Phi)T + Qk */
      matadd(matmul(matmul(state->Phi,state->P,interm1),mattrn(state->Phi,interm0),state->M),Qk,state->M);
      /* Kk = Mk * (H)T * (H * Mk * (H)T + Rk)^-1 */
      mattrn(state->H,interm2);
      matsca(matmul(state->M,interm2,state->K),state->K,
             1.0 / (matmul(matmul(state->H,state->M,interm3),interm2,interm4)->vals[0] + state->R));
      /* Pk = (I - Kk * H) * Mk) */
      matmul(matsub(matidn(interm0),matmul(state->K,state->H,interm1),interm0),state->M,state->P);

      /* now update our estimate of the state */
      matdup(state->A,interm2);
      matadd(matmul(state->Phi,interm2,state->A),matsca(state->K,interm5,measurement 
             - matmul(matmul(state->H,state->Phi,interm3),interm2,interm4)->vals[0]),state->A);

      free(interm0);
      free(interm1);
      free(interm2);
      free(interm3);
      free(interm4);
      free(interm5);
   }
   free(Qk);
}

void init_second_order_polynomial_kalman_filter(sopkfs * state) {
   matzro(state->Phi = DALLOCMAT(3,3));
   matzro(state->H = DALLOCMAT(1,3));
   matzro(state->P = DALLOCMAT(3,3));
   matzro(state->M = DALLOCMAT(3,3));
   matzro(state->K = DALLOCMAT(3,1));
   matzro(state->A = DALLOCMAT(3,1));
}

void destroy_second_order_polynomial_kalman_filter(sopkfs * state) {
   free(state->Phi);
   free(state->H);
   free(state->P);
   free(state->M);
   free(state->K);
   free(state->A);
}

void step_second_order_polynomial_kalman_filter(double measurement, double elapsed, sopkfs * state) {

   /* prepare the matrix Qk for this update */
   mat Qk = DALLOCMAT(3,3);
   Qk->vals[0] = (elapsed*elapsed*elapsed*elapsed*elapsed) / 20.0;
   Qk->vals[1] = Qk->vals[3] = (elapsed*elapsed*elapsed*elapsed) / 8.0;
   Qk->vals[2] = Qk->vals[4] = Qk->vals[6] = (elapsed*elapsed*elapsed) / 6.0;
   Qk->vals[5] = Qk->vals[7] = (elapsed*elapsed) / 2.0;
   Qk->vals[8] = elapsed;
   matsca(Qk,Qk,state->Q);

   {
      mat interm0 = DALLOCMAT(3,3),
          interm1 = DALLOCMAT(3,3),
          interm2 = DALLOCMAT(3,1),
          interm3 = DALLOCMAT(1,3),
          interm4 = DALLOCMAT(1,1),
          interm5 = DALLOCMAT(3,1);

      /* perform the Riccati equations to update the kalman gain */

      /* Mk = Phi * Pk-1 * (Phi)T + Qk */
      matadd(matmul(matmul(state->Phi,state->P,interm1),mattrn(state->Phi,interm0),state->M),Qk,state->M);
      /* Kk = Mk * (H)T * (H * Mk * (H)T + Rk)^-1 */
      mattrn(state->H,interm2);
      matsca(matmul(state->M,interm2,state->K),state->K,
             1.0 / (matmul(matmul(state->H,state->M,interm3),interm2,interm4)->vals[0] + state->R));
      /* Pk = (I - Kk * H) * Mk) */
      matmul(matsub(matidn(interm0),matmul(state->K,state->H,interm1),interm0),state->M,state->P);

      /* now update our estimate of the state */
      matdup(state->A,interm2);
      matadd(matmul(state->Phi,interm2,state->A),matsca(state->K,interm5,measurement 
             - matmul(matmul(state->H,state->Phi,interm3),interm2,interm4)->vals[0]),state->A);

      free(interm0);
      free(interm1);
      free(interm2);
      free(interm3);
      free(interm4);
      free(interm5);
   }
   free(Qk);

}

void step_gravity_informed_first_order_polynomial_kalman_filter(double measurement, double elapsed, fopkfs * state) {

   /* prepare the matrix Qk for this update */
   mat Qk = DALLOCMAT(2,2);
   Qk->vals[0] = (elapsed*elapsed*elapsed) / 3.0;
   Qk->vals[1] = Qk->vals[2] = (elapsed*elapsed) / 2.0;
   Qk->vals[3] = elapsed;
   matsca(Qk,Qk,state->Q);

   {
      mat interm0 = DALLOCMAT(2,2),
          interm1 = DALLOCMAT(2,2),
          interm2 = DALLOCMAT(2,1),
          interm3 = DALLOCMAT(1,2),
          interm4 = DALLOCMAT(1,1),
          interm5 = DALLOCMAT(2,1);

      /* perform the Riccati equations to update the kalman gain */

      /* Mk = Phi * Pk-1 * (Phi)T + Qk */
      matadd(matmul(matmul(state->Phi,state->P,interm1),mattrn(state->Phi,interm0),state->M),Qk,state->M);
      /* Kk = Mk * (H)T * (H * Mk * (H)T + Rk)^-1 */
      mattrn(state->H,interm2);
      matsca(matmul(state->M,interm2,state->K),state->K,
             1.0 / (matmul(matmul(state->H,state->M,interm3),interm2,interm4)->vals[0] + state->R));
      /* Pk = (I - Kk * H) * Mk) */
      matmul(matsub(matidn(interm0),matmul(state->K,state->H,interm1),interm0),state->M,state->P);

      /* now update our estimate of the state */
      matdup(state->A,interm2);
      matadd(matmul(state->Phi,interm2,state->A),matsca(state->K,interm5,measurement  /* V gravity term V */
             - matmul(matmul(state->H,state->Phi,interm3),interm2,interm4)->vals[0] + 16.1*elapsed*elapsed),state->A);

      /* add other gravity terms */
      state->A->vals[0] -= 16.1*elapsed*elapsed;
      state->A->vals[1] -= 32.2*elapsed;

      free(interm0);
      free(interm1);
      free(interm2);
      free(interm3);
      free(interm4);
      free(interm5);
   }
   free(Qk);
}
