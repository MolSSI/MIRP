#pragma once

#include <gmp.h>
#include <mpfr.h>

#ifdef __cplusplus
extern "C" {
#endif

void mirp_single_eri(double * result,
                       int l1, int m1, int n1, double alpha1, double A[3],
                       int l2, int m2, int n2, double alpha2, double B[3],
                       int l3, int m3, int n3, double alpha3, double C[3],
                       int l4, int m4, int n4, double alpha4, double D[3]);

void mirp_single_eri_mp(mpfr_t result,
                          int l1, int m1, int n1, mpfr_t alpha1, mpfr_t A[3],
                          int l2, int m2, int n2, mpfr_t alpha2, mpfr_t B[3],
                          int l3, int m3, int n3, mpfr_t alpha3, mpfr_t C[3],
                          int l4, int m4, int n4, mpfr_t alpha4, mpfr_t D[3],
                          mpfr_prec_t working_precision);

void mirp_single_eri_interval(arb_t result,
                              int l1, int m1, int n1, arb_t alpha1, arb_t A[3],
                              int l2, int m2, int n2, arb_t alpha2, arb_t B[3],
                              int l3, int m3, int n3, arb_t alpha3, arb_t C[3],
                              int l4, int m4, int n4, arb_t alpha4, arb_t D[3],
                              slong working_prec);

#ifdef __cplusplus
}
#endif

