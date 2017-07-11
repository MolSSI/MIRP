/*! \file
 *
 * \brief Calculation of the Gaussian Product Theorem terms
 */

#pragma once

#include <gmp.h>
#include <mpfr.h>
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
 * to be arrays of length three.
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


/*! \brief Computes factors for the Gaussian Product Theorem using arbitrary precision
 *
 *  \copydetails mirp_gpt_double
 *  \param [in] working_prec The working precision to use in the calculation
 */
void mirp_gpt_mp(const mpfr_t alpha1, const mpfr_t alpha2,
                 const mpfr_t * A, const mpfr_t * B,
                 mpfr_t gamma, mpfr_t * P,
                 mpfr_t * PA, mpfr_t * PB,
                 mpfr_t AB2,
                 mpfr_prec_t working_prec);


/*! \brief Computes factors for the Gaussian Product Theorem using interval arithmetic
 *
 *  \copydetails mirp_gpt_double
 *  \param [in] working_prec The working precision to use in the calculation
 */
void mirp_gpt_interval(const arb_t alpha1, const arb_t alpha2,
                       const arb_t * A, const arb_t * B,
                       arb_t gamma,
                       arb_t * P, arb_t * PA, arb_t * PB,
                       arb_t AB2,
                       slong working_prec);

#ifdef __cplusplus
}
#endif

