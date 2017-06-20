/*! \file
 *
 * \brief Calculation of the boys function in double precision
 *        and arbitrary precision
 */

#pragma once

#include <gmp.h>
#include <mpfr.h>
#include <arb.h>

#ifdef __cplusplus
extern "C" {
#endif

void mirp_boys_double(double *F, int m, double t);

void mirp_boys_mp(mpfr_t *F, int m, const mpfr_t t, mpfr_prec_t working_prec);

void mirp_boys_interval(arb_t *F, int m, const arb_t t, slong working_prec);


#ifdef __cplusplus
}
#endif

