/*! \file
 *
 * \brief Calculation of the boys function
 */

#pragma once

#include <arb.h>

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Computes the Boys function using double precision
 *
 * See \ref boys_function
 *
 * \warning \p F must be large enough to hold (\p m + 1) values, since
 *             this is computing from zero to m.
 *
 * \param [out] F The computed values of the Boys function
 * \param [in] m The maximum order to calculate
 * \param [in] t The value at which to evaluate
 */
void mirp_boys_d(double *F, int m, double t);


/*! \brief Computes the Boys function using interval arithmetic
 *
 * \copydetails mirp_boys_d
 * \param [in] working_prec The working precision (binary digits/bits)
 *             to use in the calculation
 */
void mirp_boys(arb_ptr F, int m, const arb_t t, slong working_prec);


/*! \brief Compute the Boys function to a target precision
 *
 * This function seeks to calculate the Boys function to a given
 * precision. This works by increasing the working precision
 * until the result is the desired precision (\p target_prec).
 *
 * If the inputs are not of high enough precision to reach the target
 * precision, the function will return nonzero, indicating that
 * reaching the target precision is impossible. Otherwise, the
 * function returns zero.
 *
 * Note that the result may be of higher precision that the
 * desired target precision.
 *
 * \param [out] F The computed values of the Boys function
 * \param [in] m The maximum order to calculate
 * \param [in] t The value at which to evaluate
 * \param [in] target_prec The desired precision we want in the computed values
 *
 * \return 0 If the function is successful, nonzero if the target
 *           precision can not be reached.
 */
int mirp_boys_target(arb_ptr F, int m, const arb_t t, slong target_prec);


/*! \brief Compute the Boys function to a target precision (string input)
 *
 * This function seeks to calculate the Boys function to a given
 * precision. This works by increasing the working precision
 * until the result is the desired precision (\p target_prec).
 *
 * The input is interpreted to be to infinite precision (ie,
 * there are an infinite number of trailing zeroes after the
 * decimal point).
 *
 * Note that the result may be of higher precision that the
 * desired target precision.
 *
 * \param [out] F The computed values of the Boys function
 * \param [in] m The maximum order to calculate
 * \param [in] t The value at which to evaluate
 * \param [in] target_prec The desired precision we want in the computed values
 */
void mirp_boys_target_str(arb_ptr F, int m, const char * t, slong target_prec);


/*! \brief Computes the Boys function to double precision using interval arithmetic
 *
 * This function takes double precision as input and returns double precision
 * as output. Internally, it uses interval arithmetic to ensure that no
 * precision is lost
 *
 * \param [out] F The computed values of the Boys function
 * \param [in]  m The maximum order to calculate
 * \param [in]  t Teh value at which to evaluate
 */
void mirp_boys_exact(double *F, int m, double t);


#ifdef __cplusplus
}
#endif

