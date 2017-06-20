/*! \file
 *
 * \brief Calculation of the Gaussian Product Theorem terms
 *        in double precision and arbitrary precision
 */

#pragma once

#include <gmp.h>
#include <mpfr.h>

#ifdef __cplusplus
extern "C" {
#endif

void mirp_gpt(double alpha1, double alpha2, double A[3], double B[3],
                double * gamma, double P[3], double PA[3], double PB[3],
                double * AB2);


void mirp_gpt_mp(mpfr_t alpha1, mpfr_t alpha2,
                   mpfr_t A[3], mpfr_t B[3],
                   mpfr_t gamma,
                   mpfr_t P[3], mpfr_t PA[3], mpfr_t PB[3],
                   mpfr_t AB2,
                   mpfr_prec_t working_prec);

#ifdef __cplusplus
}
#endif

