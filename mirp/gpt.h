/*! \file
 *
 * \brief Calculation of the Gaussian Product Theorem terms
 */

#pragma once

#include <arb.h>

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Computes factors for the Gaussian Product Theorem using double precision
 * 
 * The input parameters \p A and \p B, are expected to be arrays
 * of 3 elements.
 *
 * The output parameters \p P, \p PA, and \p PB are also expected
 * to be arrays 3 elements.
 *
 * See \ref gaussian_product_theorem
 *
 * \param [in]  alpha1  Exponent of the first gaussian
 * \param [in]  alpha2  Exponent of the second gaussian
 * \param [in]  A       XYZ coordinates of the first gaussian.
 * \param [in]  B       XYZ coordinates of the second gaussian
 * \param [out] gamma   Combined exponent (\p alpha1 + \p alpha2)
 * \param [out] P       Coordinates of the new gaussian
 * \param [out] PA      XYZ distances between the new gaussian and the first gaussian
 * \param [out] PB      XYZ distances between the new gaussian and the second gaussian
 * \param [out] AB2     Total distances between the first and second gaussians
 */
void mirp_gpt_double(double alpha1, double alpha2,
                     const double * A, const double * B,
                     double * gamma, double * P,
                     double * PA, double * PB,
                     double * AB2);


/*! \brief Computes factors for the Gaussian Product Theorem using interval arithmetic
 *
 * \copydetails mirp_gpt_double
 * \param [in] working_prec The working precision (binary digits/bits) to use in the calculation
 */
void mirp_gpt_interval(const arb_t alpha1, const arb_t alpha2,
                       arb_srcptr A, arb_srcptr B,
                       arb_t gamma,
                       arb_ptr P, arb_ptr PA, arb_ptr PB,
                       arb_t AB2,
                       slong working_prec);

#ifdef __cplusplus
}
#endif

