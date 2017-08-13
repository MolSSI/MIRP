/*! \file
 *
 * \brief Calculation of the Gaussian Product Theorem terms
 */

#include "mirp/gpt.h"

void mirp_gpt_d(double alpha1, double alpha2,
                     const double * A, const double * B,
                     double * gamma, double * P,
                     double * PA, double * PB,
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


void mirp_gpt(const arb_t alpha1, const arb_t alpha2,
                       arb_srcptr A, arb_srcptr B,
                       arb_t gamma, arb_ptr P,
                       arb_ptr PA, arb_ptr PB,
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

    arb_mul(tmp1, alpha1, A+0,   working_prec);
    arb_mul(tmp2, alpha2, B+0,   working_prec);
    arb_add(P+0, tmp1,   tmp2,   working_prec);
    arb_div(P+0, P+0,   gamma,  working_prec);

    arb_mul(tmp1, alpha1, A+1,   working_prec);
    arb_mul(tmp2, alpha2, B+1,   working_prec);
    arb_add(P+1, tmp1,   tmp2,   working_prec);
    arb_div(P+1, P+1,   gamma,  working_prec);

    arb_mul(tmp1, alpha1, A+2,   working_prec);
    arb_mul(tmp2, alpha2, B+2,   working_prec);
    arb_add(P+2, tmp1,   tmp2,   working_prec);
    arb_div(P+2, P+2,   gamma,  working_prec);

    /*
     *
     * PA[0] = P[0] - A[0], etc
     * PB[0] = P[0] - B[0], etc
     */
    arb_sub(PA+0, P+0, A+0, working_prec);
    arb_sub(PA+1, P+1, A+1, working_prec);
    arb_sub(PA+2, P+2, A+2, working_prec);
    arb_sub(PB+0, P+0, B+0, working_prec);
    arb_sub(PB+1, P+1, B+1, working_prec);
    arb_sub(PB+2, P+2, B+2, working_prec);

    /*
     * AB2 = (A[0]-B[0])*(A[0]-B[0]) + (A[1]-B[1])*(A[1]-B[1]) + (A[2]-B[2])*(A[2]-B[2]);
     */
    arb_sub(tmp1, A+0, B+0, working_prec);
    arb_sub(tmp2, A+1, B+1, working_prec);
    arb_sub(tmp3, A+2, B+2, working_prec);
    arb_mul(AB2,  tmp1, tmp1, working_prec);
    arb_addmul(AB2, tmp2, tmp2, working_prec);
    arb_addmul(AB2, tmp3, tmp3, working_prec);

    /* clean up temporaries */
    arb_clear(tmp1);
    arb_clear(tmp2);
    arb_clear(tmp3);
}
