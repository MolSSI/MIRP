#pragma once

#include <gmp.h>
#include <mpfr.h>
#include <arb.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

void mirp_single_eri_double(double * result,
                            const int * lmn1, const double * A, double alpha1,
                            const int * lmn2, const double * B, double alpha2,
                            const int * lmn3, const double * C, double alpha3,
                            const int * lmn4, const double * D, double alpha4);

size_t mirp_prim_eri_double(double * result,
                            int am1, const double * A, double alpha1,
                            int am2, const double * B, double alpha2,
                            int am3, const double * C, double alpha3,
                            int am4, const double * D, double alpha4);

size_t mirp_eri_double(double * result,
                       int am1, const double * A, int nprim1, int ngeneral1, const double * alpha1, const double * coeff1,
                       int am2, const double * B, int nprim2, int ngeneral2, const double * alpha2, const double * coeff2,
                       int am3, const double * C, int nprim3, int ngeneral3, const double * alpha3, const double * coeff3,
                       int am4, const double * D, int nprim4, int ngeneral4, const double * alpha4, const double * coeff4);

void mirp_single_eri_mp(mpfr_t result,
                        const int * lmn1, const mpfr_t * A, const mpfr_t alpha1,
                        const int * lmn2, const mpfr_t * B, const mpfr_t alpha2,
                        const int * lmn3, const mpfr_t * C, const mpfr_t alpha3,
                        const int * lmn4, const mpfr_t * D, const mpfr_t alpha4,
                        mpfr_prec_t working_precision);

size_t mirp_prim_eri_mp(mpfr_t * result,
                        int am1, const mpfr_t alpha1, const mpfr_t * A,
                        int am2, const mpfr_t alpha2, const mpfr_t * B,
                        int am3, const mpfr_t alpha3, const mpfr_t * C,
                        int am4, const mpfr_t alpha4, const mpfr_t * D);

size_t mirp_eri_mp(mpfr_t * result,
                   int am1, const mpfr_t * A, int nprim1, int ngeneral1, const mpfr_t * alpha1, const mpfr_t * coeff1,
                   int am2, const mpfr_t * B, int nprim2, int ngeneral2, const mpfr_t * alpha2, const mpfr_t * coeff2,
                   int am3, const mpfr_t * C, int nprim3, int ngeneral3, const mpfr_t * alpha3, const mpfr_t * coeff3,
                   int am4, const mpfr_t * D, int nprim4, int ngeneral4, const mpfr_t * alpha4, const mpfr_t * coeff4,
                   mpfr_prec_t working_prec);

void mirp_single_eri_interval(arb_t result,
                              const int * lmn1, const arb_t alpha1, const arb_t * A,
                              const int * lmn2, const arb_t alpha2, const arb_t * B,
                              const int * lmn3, const arb_t alpha3, const arb_t * C,
                              const int * lmn4, const arb_t alpha4, const arb_t * D,
                              slong working_prec);

size_t mirp_prim_eri_interval(arb_t * result,
                        int am1, const arb_t * A, const arb_t alpha1,
                        int am2, const arb_t * B, const arb_t alpha2,
                        int am3, const arb_t * C, const arb_t alpha3,
                        int am4, const arb_t * D, const arb_t alpha4,
                        slong working_prec);

size_t mirp_eri_interval(arb_t * result,
                   int am1, const arb_t * A, int nprim1, int ngeneral1, const arb_t * alpha1, const arb_t * coeff1,
                   int am2, const arb_t * B, int nprim2, int ngeneral2, const arb_t * alpha2, const arb_t * coeff2,
                   int am3, const arb_t * C, int nprim3, int ngeneral3, const arb_t * alpha3, const arb_t * coeff3,
                   int am4, const arb_t * D, int nprim4, int ngeneral4, const arb_t * alpha4, const arb_t * coeff4,
                   slong working_prec);

#ifdef __cplusplus
}
#endif

