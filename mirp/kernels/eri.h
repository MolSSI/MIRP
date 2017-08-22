/*! \file
 *
 * \brief Calculation of electron repulsion integrals
 */

#pragma once

#include <arb.h>
#include "mirp/kernels/integral4_wrappers.h"

#ifdef __cplusplus
extern "C" {
#endif


/*! \brief Computes a single cartesian electron repulsion integral
 *         (double precision)
 *
 * The \p lmn paramters are the exponents on x, y, and z, with the total
 * angular momentum being the sum of the three components. For example,
 * { 1, 2, 0 } is the \f$f_{xy^2}\f$ cartesian integral.
 *
 * \param [out] integral
 *              Resulting integral integral
 * \param [in]  lmn1,lmn2,lmn3,lmn4
 *              Exponents of x, y, and z that signify angular momentum. Required
 *              to be 3 elements.
 * \param [in]  A,B,C,D
 *              XYZ coordinates of the four centers (each of length 3)
 * \param [in]  alpha1,alpha2,alpha3,alpha4
 *              Exponents of the gaussian on the four centers
 */
void mirp_eri_single_d(double * integral,
                       const int * lmn1, const double * A, double alpha1,
                       const int * lmn2, const double * B, double alpha2,
                       const int * lmn3, const double * C, double alpha3,
                       const int * lmn4, const double * D, double alpha4);


/*! \brief Computes a single cartesian electron repulsion integral
 *         (interval arithmetic)
 *
 * \copydetails mirp_eri_single_d
 * \param [in] working_prec The working precision (binary digits/bits) to use in the calculation
 */
void mirp_eri_single(arb_t integral,
                     const int * lmn1, arb_srcptr A, const arb_t alpha1,
                     const int * lmn2, arb_srcptr B, const arb_t alpha2,
                     const int * lmn3, arb_srcptr C, const arb_t alpha3,
                     const int * lmn4, arb_srcptr D, const arb_t alpha4,
                     slong working_prec);


/****************************************
 * Wrappings
 ****************************************/
MIRP_WRAP_SHELL4(eri)
MIRP_WRAP_SHELL4_D(eri)
MIRP_WRAP_SINGLE4_TARGET(eri)
MIRP_WRAP_SINGLE4_TARGET_STR(eri)
MIRP_WRAP_SINGLE4_EXACT(eri)

MIRP_WRAP_TARGET(eri)
MIRP_WRAP_EXACT(eri)


#ifdef __cplusplus
}
#endif

