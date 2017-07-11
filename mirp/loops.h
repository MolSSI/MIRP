/*! \file
 *
 * \brief Helpers for looping over cartesian functions and primitives
 */

#pragma once

#include <gmp.h>
#include <mpfr.h>
#include <arb.h>


#ifdef __cplusplus
extern "C" {
#endif

/*! \brief A callback for a function that takes 4 gaussian basis functions */
typedef void (*cb_4gaussians_double)(double *,
                                     const int *, const double *, double,
                                     const int *, const double *, double,
                                     const int *, const double *, double,
                                     const int *, const double *, double);

typedef void (*cb_4shells_double)(double *,
                                  int, const double *, double,
                                  int, const double *, double,
                                  int, const double *, double,
                                  int, const double *, double);

typedef void (*cb_4gaussians_mp)(mpfr_t,
                                 const int *, const mpfr_t *, const mpfr_t,
                                 const int *, const mpfr_t *, const mpfr_t,
                                 const int *, const mpfr_t *, const mpfr_t,
                                 const int *, const mpfr_t *, const mpfr_t,
                                 mpfr_prec_t);

typedef void (*cb_4shells_mp)(mpfr_t *,
                              int, const mpfr_t *, const mpfr_t,
                              int, const mpfr_t *, const mpfr_t,
                              int, const mpfr_t *, const mpfr_t,
                              int, const mpfr_t *, const mpfr_t,
                              mpfr_prec_t);

typedef void (*cb_4gaussians_interval)(arb_t,
                                       const int *, const arb_t *, const arb_t,
                                       const int *, const arb_t *, const arb_t,
                                       const int *, const arb_t *, const arb_t,
                                       const int *, const arb_t *, const arb_t,
                                       slong);

typedef void (*cb_4shells_interval)(arb_t *,
                                    int, const arb_t *, const arb_t,
                                    int, const arb_t *, const arb_t,
                                    int, const arb_t *, const arb_t,
                                    int, const arb_t *, const arb_t,
                                    slong);



void mirp_prim_cartloop4_double(double * result,
                                int am1, const double * A, double alpha1,
                                int am2, const double * B, double alpha2,
                                int am3, const double * C, double alpha3,
                                int am4, const double * D, double alpha4,
                                cb_4gaussians_double cb);

void mirp_cartloop4_double(double * result,
                           int am1, const double * A, int nprim1, int ngeneral1, const double * alpha1, const double * coeff1,
                           int am2, const double * B, int nprim2, int ngeneral2, const double * alpha2, const double * coeff2,
                           int am3, const double * C, int nprim3, int ngeneral3, const double * alpha3, const double * coeff3,
                           int am4, const double * D, int nprim4, int ngeneral4, const double * alpha4, const double * coeff4,
                           cb_4shells_double cb);



void mirp_prim_cartloop4_mp(mpfr_t * result,
                            int am1, const mpfr_t * A, const mpfr_t alpha1,
                            int am2, const mpfr_t * B, const mpfr_t alpha2,
                            int am3, const mpfr_t * C, const mpfr_t alpha3,
                            int am4, const mpfr_t * D, const mpfr_t alpha4,
                            mpfr_prec_t working_prec, cb_4gaussians_mp cb);

void mirp_cartloop4_mp(mpfr_t * result,
                       int am1, const mpfr_t * A, int nprim1, int ngeneral1, const mpfr_t * alpha1, const mpfr_t * coeff1,
                       int am2, const mpfr_t * B, int nprim2, int ngeneral2, const mpfr_t * alpha2, const mpfr_t * coeff2,
                       int am3, const mpfr_t * C, int nprim3, int ngeneral3, const mpfr_t * alpha3, const mpfr_t * coeff3,
                       int am4, const mpfr_t * D, int nprim4, int ngeneral4, const mpfr_t * alpha4, const mpfr_t * coeff4,
                       mpfr_prec_t working_prec, cb_4shells_mp cb);


void mirp_prim_cartloop4_interval(arb_t * result,
                                  int am1, const arb_t * A, const arb_t alpha1,
                                  int am2, const arb_t * B, const arb_t alpha2,
                                  int am3, const arb_t * C, const arb_t alpha3,
                                  int am4, const arb_t * D, const arb_t alpha4,
                                  slong working_prec, cb_4gaussians_interval cb);

void mirp_cartloop4_interval(arb_t * result,
                             int am1, const arb_t * A, int nprim1, int ngeneral1, const arb_t * alpha1, const arb_t * coeff1,
                             int am2, const arb_t * B, int nprim2, int ngeneral2, const arb_t * alpha2, const arb_t * coeff2,
                             int am3, const arb_t * C, int nprim3, int ngeneral3, const arb_t * alpha3, const arb_t * coeff3,
                             int am4, const arb_t * D, int nprim4, int ngeneral4, const arb_t * alpha4, const arb_t * coeff4,
                             slong working_prec, cb_4shells_interval cb);


#ifdef __cplusplus
}
#endif

