
/* implimentations of the filters presented in chapter 4. We are kind of at
 * the point where the filter needs to be constructed mathmatically for
 * the specific task we are approaching. This is interesting, as I would 
 * have assumed that the system was more general. It is general in the sense
 * that the Riccati equations do not change, but at the same time the matrix
 * F and matrix G seem to need to be reengineered for the given problem.
 * This makes sense given that the Kalman filter is basically the same as the
 * least-squares filter without G and Q. It is preferred because it can be
 * customized and augmented with additional information.
 */

#ifndef __TANG_POLYNOMIAL_KALMAN_FILTERING__

#define __TANG_POLYNOMIAL_KALMAN_FILTERING__
#include "numerical_basics.h"

/* holds state information for the filter. All values are 
 * numbers instead of matrices because if they were matrices
 * they would be matrices with a single dimension.
 *
 * Phi - fundamental update matrix
 * H   - measurement update matrix
 * R   - measurement noise matrix
 * Q   - process noise spectral density
 * P   - covariance of state error
 * M   - covariance of state error before update
 *
 * a0  - position estimate
 */
typedef struct zeroth_order_polynomial_kalman_filter_state {
   double a0;
   double Phi;
   double H;
   double R;
   double Q;
   double P;
   double M;
   double K;
} zopkfs;

/* update the kalman filter with the given measurement and timestep. Unlike the
 * filter designs in the previous files, elapsed here is the time distance 
 * between the last measurement provided to the filter and the current one.
 * Ex: if you receive a measurement every x seconds, x would be the value
 *     you pass to elapsed every time you recieve a measurement.
 */
void step_zeroth_order_polynomial_kalman_filter(double measurement, double elapsed, zopkfs * state);

/* holds state information for the filter.
 *
 * Phi - fundamental update matrix [2x2]
 * H   - measurement update matrix [1x2]
 * R   - measurement noise matrix  [1x1]
 * Q   - process noise spectral density 
 * P   - covariance of state error [2x2]
 * M   - covariance of prior error [2x2]
 * K   - kalman gain               [2x1]
 *
 * A[0]  - position estimate
 * A[1]  - velocity estimate
 */
typedef struct first_order_polynomial_kalman_filter_state {
   mat A;
   mat Phi;
   mat H;
   mat P;
   mat M;
   mat K;
   double R;
   double Q;
} fopkfs;

/* allocate and zero the matrices to be used by the filter */
void init_first_order_polynomial_kalman_filter(fopkfs * state);

/* deallocate the matrices used by the filter */
void destroy_first_order_polynomial_kalman_filter(fopkfs * state);

/* update the kalman filter with the given measurement and timestep. */
void step_first_order_polynomial_kalman_filter(double measurement, double elapsed, fopkfs * state);

/* holds state information for the filter.
 *
 * Phi - fundamental update matrix [3x3]
 * H   - measurement update matrix [1x3]
 * R   - measurement noise matrix  [1x1]
 * Q   - process noise spectral density 
 * P   - covariance of state error [3x3]
 * M   - covariance of prior error [3x3]
 * K   - kalman gain               [2x1]
 *
 * A[0]  - position estimate
 * A[1]  - velocity estimate
 * A[2]  - acceleration estimate
 */
typedef struct second_order_polynomial_kalman_filter_state {
   mat A;
   mat Phi;
   mat H;
   mat P;
   mat M;
   mat K;
   double R;
   double Q;
} sopkfs;

/* allocate and zero the matrices to be used by the filter */
void init_second_order_polynomial_kalman_filter(sopkfs * state);

/* deallocate the matrices used by the filter */
void destroy_second_order_polynomial_kalman_filter(sopkfs * state);

/* update the kalman filter with the given measurement and timestep. */
void step_second_order_polynomial_kalman_filter(double measurement, double elapsed, sopkfs * state);

/* a first order kalman filter that takes into account the influence
 * of gravity. The rest of the math and auxillary variables do not
 * change, so just the update function needs to be changed. In reality
 * it would be better to add a G tern to the struct and incorporate
 * it into the full equation, but I will save that for the future... */
void step_gravity_informed_first_order_polynomial_kalman_filter(double measurement, double elapsed, fopkfs * state);





#endif
