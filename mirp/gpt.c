/*! \file
 *
 * \brief Calculation of the Gaussian Product Theorem terms
 *        in double precision and arbitrary precision
 */

#include "mirp/gpt.h"

void mirp_gpt(double alpha1, double alpha2, double A[3], double B[3],
              double * gamma, double P[3], double PA[3], double PB[3],
              double * AB2)
{
    *gamma = alpha1 + alpha2;

    P[0] = (alpha1*A[0] + alpha2*B[0]) / (*gamma);
    P[1] = (alpha1*A[1] + alpha2*B[1]) / (*gamma);
    P[2] = (alpha1*A[2] + alpha2*B[2]) / (*gamma);

    PA[0] = P[0] - A[0];
    PA[1] = P[1] - A[1];
    PA[2] = P[2] - A[2];

    PB[0] = P[0] - B[0];
    PB[1] = P[1] - B[1];
    PB[2] = P[2] - B[2];

    *AB2 = (A[0]-B[0])*(A[0]-B[0])
         + (A[1]-B[1])*(A[1]-B[1])
         + (A[2]-B[2])*(A[2]-B[2]);
}


void mirp_gpt_mp(mpfr_t alpha1, mpfr_t alpha2,
                 mpfr_t A[3], mpfr_t B[3],
                 mpfr_t gamma,
                 mpfr_t P[3], mpfr_t PA[3], mpfr_t PB[3],
                 mpfr_t AB2,
                 mpfr_prec_t working_prec)
{
    /* Temporary data */
    mpfr_t tmp1, tmp2, tmp3;
    mpfr_inits2(working_prec, tmp1, tmp2, tmp3, (mpfr_ptr)0);


    /*
     * gamma = alpha1 + alpha2
     *
     * P[0] = (alpha1*A[0] + alpha2*B[0]) / gamma
     * and similar for P[1], P[2]
     */
    mpfr_add(gamma, alpha1, alpha2, MPFR_RNDN);

    mpfr_mul(tmp1, alpha1, A[0],   MPFR_RNDN);
    mpfr_mul(tmp2, alpha2, B[0],   MPFR_RNDN);
    mpfr_add(P[0], tmp1,   tmp2,   MPFR_RNDN);
    mpfr_div(P[0], P[0],   gamma,  MPFR_RNDN);

    mpfr_mul(tmp1, alpha1, A[1],   MPFR_RNDN);
    mpfr_mul(tmp2, alpha2, B[1],   MPFR_RNDN);
    mpfr_add(P[1], tmp1,   tmp2,   MPFR_RNDN);
    mpfr_div(P[1], P[1],   gamma,  MPFR_RNDN);

    mpfr_mul(tmp1, alpha1, A[2],   MPFR_RNDN);
    mpfr_mul(tmp2, alpha2, B[2],   MPFR_RNDN);
    mpfr_add(P[2], tmp1,   tmp2,   MPFR_RNDN);
    mpfr_div(P[2], P[2],   gamma,  MPFR_RNDN);

    /*
     *
     * PA[0] = P[0] - A[0], etc
     * PB[0] = P[0] - B[0], etc
     */
    mpfr_sub(PA[0], P[0], A[0], MPFR_RNDN);
    mpfr_sub(PA[1], P[1], A[1], MPFR_RNDN);
    mpfr_sub(PA[2], P[2], A[2], MPFR_RNDN);
    mpfr_sub(PB[0], P[0], B[0], MPFR_RNDN);
    mpfr_sub(PB[1], P[1], B[1], MPFR_RNDN);
    mpfr_sub(PB[2], P[2], B[2], MPFR_RNDN);

    /*
     * AB2 = (A[0]-B[0])*(A[0]-B[0]) + (A[1]-B[1])*(A[1]-B[1]) + (A[2]-B[2])*(A[2]-B[2]);
     */
    mpfr_sub(tmp1, A[0], B[0], MPFR_RNDN);
    mpfr_sub(tmp2, A[1], B[1], MPFR_RNDN);
    mpfr_sub(tmp3, A[2], B[2], MPFR_RNDN);
    mpfr_mul(AB2,  tmp1, tmp1, MPFR_RNDN);
    mpfr_fma(AB2,  tmp2, tmp2, AB2, MPFR_RNDN);
    mpfr_fma(AB2,  tmp3, tmp3, AB2, MPFR_RNDN);

    /* clean up temporaries */
    mpfr_clears(tmp1, tmp2, tmp3, (mpfr_ptr)0);
}

void mirp_gpt_interval(arb_t alpha1, arb_t alpha2,
                       arb_t A[3], arb_t B[3],
                       arb_t gamma,
                       arb_t P[3], arb_t PA[3], arb_t PB[3],
                       arb_t AB2,
                       slong working_prec)
{
    /* Temporary data */
    arb_t tmp1, tmp2, tmp3;
    arb_init(tmp1);
    arb_init(tmp2);
    arb_init(tmp3);

    /*
     * gamma = alpha1 + alpha2
     *
     * P[0] = (alpha1*A[0] + alpha2*B[0]) / gamma
     * and similar for P[1], P[2]
     */
    arb_add(gamma, alpha1, alpha2, working_prec);

    arb_mul(tmp1, alpha1, A[0],   working_prec);
    arb_mul(tmp2, alpha2, B[0],   working_prec);
    arb_add(P[0], tmp1,   tmp2,   working_prec);
    arb_div(P[0], P[0],   gamma,  working_prec);

    arb_mul(tmp1, alpha1, A[1],   working_prec);
    arb_mul(tmp2, alpha2, B[1],   working_prec);
    arb_add(P[1], tmp1,   tmp2,   working_prec);
    arb_div(P[1], P[1],   gamma,  working_prec);

    arb_mul(tmp1, alpha1, A[2],   working_prec);
    arb_mul(tmp2, alpha2, B[2],   working_prec);
    arb_add(P[2], tmp1,   tmp2,   working_prec);
    arb_div(P[2], P[2],   gamma,  working_prec);

    /*
     *
     * PA[0] = P[0] - A[0], etc
     * PB[0] = P[0] - B[0], etc
     */
    arb_sub(PA[0], P[0], A[0], working_prec);
    arb_sub(PA[1], P[1], A[1], working_prec);
    arb_sub(PA[2], P[2], A[2], working_prec);
    arb_sub(PB[0], P[0], B[0], working_prec);
    arb_sub(PB[1], P[1], B[1], working_prec);
    arb_sub(PB[2], P[2], B[2], working_prec);

    /*
     * AB2 = (A[0]-B[0])*(A[0]-B[0]) + (A[1]-B[1])*(A[1]-B[1]) + (A[2]-B[2])*(A[2]-B[2]);
     */
    arb_sub(tmp1, A[0], B[0], working_prec);
    arb_sub(tmp2, A[1], B[1], working_prec);
    arb_sub(tmp3, A[2], B[2], working_prec);
    arb_mul(AB2,  tmp1, tmp1, working_prec);
    arb_addmul(AB2, tmp2, tmp2, working_prec);
    arb_addmul(AB2, tmp3, tmp3, working_prec);

    /* clean up temporaries */
    arb_clear(tmp1);
    arb_clear(tmp2);
    arb_clear(tmp3);
}
