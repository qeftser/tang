
#include "kalman_filters_in_a_nonpolynomial_world.h"


void init_full_first_order_polynomial_kalman_filter(ffopkfs * state) {
   matzro(state->A   = DALLOCMAT(2,1));
   matzro(state->Phi = DALLOCMAT(2,2));
   matzro(state->H   = DALLOCMAT(1,2));
   matzro(state->P   = DALLOCMAT(2,2));
   matzro(state->M   = DALLOCMAT(2,2));
   matzro(state->K   = DALLOCMAT(2,1));
   matzro(state->G   = DALLOCMAT(2,1));
   state->R = 0.0;
   state->Q = 0.0;
}

void destroy_full_first_order_polynomial_kalman_filter(ffopkfs * state) {
   free(state->A);
   free(state->Phi);
   free(state->H);
   free(state->P);
   free(state->M);
   free(state->K);
   free(state->G);
}

void step_full_first_order_polynomial_kalman_filter(double z, double u, double elapsed, ffopkfs * state) {
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
          interm5 = DALLOCMAT(2,1),
          interm6 = DALLOCMAT(2,1),
          interm7 = DALLOCMAT(1,1);

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
      matadd(matmul(state->Phi,interm2,state->A),matsca(state->K,interm5,z 
             - matmul(matmul(state->H,state->Phi,interm3),interm2,interm4)->vals[0]
             - matmul(state->H,matsca(state->G,interm6,u),interm7)->vals[0]),state->A);

      free(interm0);
      free(interm1);
      free(interm2);
      free(interm3);
      free(interm4);
      free(interm5);
      free(interm6);
      free(interm7);
   }
   free(Qk);
}
