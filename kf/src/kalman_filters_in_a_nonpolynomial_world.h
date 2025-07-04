
/* reimplimentation of the first order polynomial discrete kalman filter that
 * includes the deterministic update term G. This is the only filter that needs
 * to be implimented to reproduce the results seen in chapter 5, so this
 * is the only filter that is implimented at this stage...                   */

#ifndef __TANG_KALMAN_FILTERS_IN_A_NONPOLYNOMIAL_WORLD__

#define __TANG_KALMAN_FILTERS_IN_A_NONPOLYNOMIAL_WORLD__
#include "numerical_basics.h"

/* holds state information for the filter.
 *
 * Phi - fundamental update matrix   [2x2]
 * H   - measurement update matrix   [1x2]
 * R   - measurement noise matrix    [1x1]
 * Q   - process noise spectral density 
 * P   - covariance of state error   [2x2]
 * M   - covariance of prior error   [2x2]
 * K   - kalman gain                 [2x1]
 * G   - deterministic update matrix [2x1]
 *
 * A[0]  - position estimate
 * A[1]  - velocity estimate
 */
typedef struct full_first_order_polynomial_kalman_filter_state {
   mat A;
   mat Phi;
   mat H;
   mat P;
   mat M;
   mat K;
   mat G;
   double R;
   double Q;
} ffopkfs;

/* allocate the matrix memory and zero all values */
void init_full_first_order_polynomial_kalman_filter(ffopkfs * state);

/* deallocate the matrix memory */
void destroy_full_first_order_polynomial_kalman_filter(ffopkfs * state);

/* perform the kalman gain update and update the estimate of
 * the filter's values.
 * z       - the measurement
 * u       - the deterministic update value
 * elapsed - time since last update
 */
void step_full_first_order_polynomial_kalman_filter(double z, double u, double elapsed, ffopkfs * state);

#endif
