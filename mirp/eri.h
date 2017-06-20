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

void mirp_single_eri_mp_str(char ** result,
                              int l1, int m1, int n1, const char * alpha1, const char * A[3],
                              int l2, int m2, int n2, const char * alpha2, const char * B[3],
                              int l3, int m3, int n3, const char * alpha3, const char * C[3],
                              int l4, int m4, int n4, const char * alpha4, const char * D[3],
                              mpfr_prec_t working_precision);

#ifdef __cplusplus
}
#endif

