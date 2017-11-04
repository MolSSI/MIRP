/*! \file
 *
 * \brief Calculation of the boys function
 */

#pragma once

#include <arb.h>

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Computes the Boys function using interval arithmetic
 *
 * See \ref boys_function
 *
 * \warning \p F must be large enough to hold (\p m + 1) values, since
 *             this is computing from zero to m.
 *
 * \param [out] F The computed values of the Boys function
 * \param [in]  m The maximum order to calculate
 * \param [in]  t The value at which to evaluate
 *
 * \param [in] working_prec The working precision (binary digits/bits)
 *                          to use in the calculation
 */
void mirp_boys(arb_ptr F, int m, const arb_t t, slong working_prec);


/*! \brief Computes the Boys function using interval arithmetic
 *         from string inputs
 *
 * \copydetails mirp_boys
 *
 * \param [in] working_prec The working precision (binary digits/bits)
 *                          to use in the calculation
 */
void mirp_boys_str(arb_ptr F, int m, const char * t, slong working_prec);


/*! \brief Computes the Boys function to exact double precision using
 *         interval arithmetic
 *
 * This function takes double precision as input and returns double precision
 * as output. Internally, it uses interval arithmetic to ensure that no
 * precision is lost
 *
 * \param [out] F The computed values of the Boys function
 * \param [in]  m The maximum order to calculate
 * \param [in]  t The value at which to evaluate
 */
void mirp_boys_exact(double *F, int m, double t);


#ifdef __cplusplus
}
#endif

