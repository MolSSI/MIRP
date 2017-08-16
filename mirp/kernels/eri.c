/*! \file
 *
 * \brief Calculation of electron repulsion integrals (interval arithmetic)
 */

#include "mirp/kernels/boys.h"
#include "mirp/kernels/eri.h"
#include "mirp/math.h"
#include "mirp/gpt.h"

static void mirp_farr(arb_ptr f,
                      int lmn1, int lmn2,
                      arb_t xyz1, arb_t xyz2,
                      slong working_prec)
{
    int i, j, k;

    arb_t tmp1, tmp2;
    arb_init(tmp1);
    arb_init(tmp2);

    _arb_vec_zero(f, lmn1 + lmn2 + 1);

    for (k = 0; k <= lmn1 + lmn2; k++)
    {
        for (i = 0; i <= MIN(k,lmn1); i++)
        {
            j = k - i;
            if (j > lmn2)
                continue;

            mirp_binomial(tmp1, lmn1, i, working_prec);
            mirp_binomial(tmp2, lmn2, j, working_prec);
            arb_mul(tmp1, tmp1, tmp2, working_prec);

            if (lmn1 - i > 0)
            {
                arb_pow_ui(tmp2, xyz1, lmn1-i, working_prec);
                arb_mul(tmp1, tmp1, tmp2, working_prec);
            }
            if (lmn2 - j > 0)
            {
                arb_pow_ui(tmp2, xyz2, lmn2-j, working_prec);
                arb_mul(tmp1, tmp1, tmp2, working_prec);
            }
            arb_add(f + k, f + k, tmp1, working_prec);
        }
    }

    arb_clear(tmp1);
    arb_clear(tmp2);
}


static void mirp_G(arb_t G, arb_t fp, arb_t fq,
                   int np, int nq, int w1, int w2,
                   arb_t gammap, arb_t gammaq, arb_t gammapq,
                   slong working_prec)
{
    arb_t tmp1, tmp2;
    arb_init(tmp1);
    arb_init(tmp2);

    arb_set_si(G, NEG1_POW(np));
    arb_mul(G, G, fp, working_prec);
    arb_mul(G, G, fq, working_prec);

    arb_fac_ui(tmp1, np, working_prec);
    arb_mul(G, G, tmp1, working_prec);
    arb_fac_ui(tmp1, nq, working_prec);
    arb_mul(G, G, tmp1, working_prec);

    arb_set_si(tmp1, w1 - np);
    arb_pow(tmp1, gammap, tmp1, working_prec);
    arb_mul(G, G, tmp1, working_prec);

    arb_set_si(tmp1, w2 - nq);
    arb_pow(tmp1, gammaq, tmp1, working_prec);
    arb_mul(G, G, tmp1, working_prec);

    arb_fac_ui(tmp1, np + nq - 2 * (w1 + w2), working_prec);
    arb_mul(G, G, tmp1, working_prec);

    arb_pow_ui(tmp1, gammapq, np + nq - 2 * (w1 + w2), working_prec);
    arb_mul(G, G, tmp1, working_prec);

    arb_fac_ui(tmp1, w1, working_prec);
    arb_fac_ui(tmp2, w2, working_prec);
    arb_mul(tmp1, tmp1, tmp2, working_prec);
    arb_fac_ui(tmp2, np - 2 * w1, working_prec);
    arb_mul(tmp1, tmp1, tmp2, working_prec);
    arb_fac_ui(tmp2, nq - 2 * w2, working_prec);
    arb_mul(tmp1, tmp1, tmp2, working_prec);
    arb_div(G, G, tmp1, working_prec);

    arb_clear(tmp1);
    arb_clear(tmp2);
}

void mirp_eri_single(arb_t integral,
                     const int * lmn1, arb_srcptr A, const arb_t alpha1,
                     const int * lmn2, arb_srcptr B, const arb_t alpha2,
                     const int * lmn3, arb_srcptr C, const arb_t alpha3,
                     const int * lmn4, arb_srcptr D, const arb_t alpha4,
                     slong working_prec)
{
    const int L_l = lmn1[0]+lmn2[0]+lmn3[0]+lmn4[0];
    const int L_m = lmn1[1]+lmn2[1]+lmn3[1]+lmn4[1];
    const int L_n = lmn1[2]+lmn2[2]+lmn3[2]+lmn4[2];
    const int L = L_l + L_m + L_n;

    arb_ptr F = _arb_vec_init(L+1);
    arb_ptr flp = _arb_vec_init(lmn1[0]+lmn2[0]+1);
    arb_ptr fmp = _arb_vec_init(lmn1[1]+lmn2[1]+1);
    arb_ptr fnp = _arb_vec_init(lmn1[2]+lmn2[2]+1);
    arb_ptr flq = _arb_vec_init(lmn3[0]+lmn4[0]+1);
    arb_ptr fmq = _arb_vec_init(lmn3[1]+lmn4[1]+1);
    arb_ptr fnq = _arb_vec_init(lmn3[2]+lmn4[2]+1);

    /* Zero the integral (we will be summing into it) */
    arb_zero(integral);

    /* Temporary variables used in constructing expressions */
    arb_t tmp1, tmp2, tmp3;
    arb_init(tmp1);
    arb_init(tmp2);
    arb_init(tmp3);


    /*************************************************
     * Calculate all the various terms from the GPT
     *************************************************/
    arb_ptr P  = _arb_vec_init(3);
    arb_ptr PA = _arb_vec_init(3);
    arb_ptr PB = _arb_vec_init(3);
    arb_ptr Q  = _arb_vec_init(3);
    arb_ptr QC = _arb_vec_init(3);
    arb_ptr QD = _arb_vec_init(3);
    arb_ptr PQ = _arb_vec_init(3);

    arb_t gammap, gammaq, gammapq, AB2, CD2, PQ2;
    arb_init(gammap);
    arb_init(gammaq);
    arb_init(gammapq);
    arb_init(AB2);
    arb_init(CD2);
    arb_init(PQ2);

    /* Gaussian Product Theorem */
    mirp_gpt(alpha1, alpha2, A, B, gammap, P, PA, PB, AB2, working_prec);
    mirp_gpt(alpha3, alpha4, C, D, gammaq, Q, QC, QD, CD2, working_prec);


    /*
     * gammapq = gammap * gammaq / (gammap + gammaq);
     * PQ[0] = P[0] - Q[0]
     * etc
     */
    arb_mul(tmp1,    gammap, gammaq, working_prec);
    arb_add(tmp2,    gammap, gammaq, working_prec);
    arb_div(gammapq, tmp1,   tmp2,   working_prec);

    arb_sub(PQ+0, P+0, Q+0, working_prec);
    arb_sub(PQ+1, P+1, Q+1, working_prec);
    arb_sub(PQ+2, P+2, Q+2, working_prec);

    /*
     * PQ2 = (P[0]-Q[0])*(P[0]-Q[0]) + (P[1]-Q[1])*(P[1]-Q[1]) + (P[2]-Q[2])*(P[2]-Q[2]);
     */
    arb_mul(PQ2, PQ+0, PQ+0, working_prec);
    arb_addmul(PQ2, PQ+1, PQ+1, working_prec);
    arb_addmul(PQ2, PQ+2, PQ+2, working_prec);


    mirp_farr(flp, lmn1[0], lmn2[0], PA+0, PB+0, working_prec);
    mirp_farr(fmp, lmn1[1], lmn2[1], PA+1, PB+1, working_prec);
    mirp_farr(fnp, lmn1[2], lmn2[2], PA+2, PB+2, working_prec);
    mirp_farr(flq, lmn3[0], lmn4[0], QC+0, QD+0, working_prec);
    mirp_farr(fmq, lmn3[1], lmn4[1], QC+1, QD+1, working_prec);
    mirp_farr(fnq, lmn3[2], lmn4[2], QC+2, QD+2, working_prec);


    /*
     *  Calculate the Boys function
     */
    arb_mul(tmp1, PQ2, gammapq, working_prec);
    mirp_boys(F, L, tmp1, working_prec);


    /*
     * G values used within the loops
     */
    arb_t Gx, Gy, Gz, Gxy, Gxyz;
    arb_init(Gx);
    arb_init(Gy);
    arb_init(Gz);
    arb_init(Gxy);
    arb_init(Gxyz);

    for(int lp = 0; lp <= lmn1[0] + lmn2[0]; lp++)
    for(int lq = 0; lq <= lmn3[0] + lmn4[0]; lq++)
    for(int u1 = 0; u1 <= (lp/2); u1++)
    for(int u2 = 0; u2 <= (lq/2); u2++)
    {
        mirp_G(Gx, flp + lp, flq + lq, lp, lq, u1, u2, gammap, gammaq, gammapq, working_prec);

        for(int mp = 0; mp <= lmn1[1] + lmn2[1]; mp++)
        for(int mq = 0; mq <= lmn3[1] + lmn4[1]; mq++)
        for(int v1 = 0; v1 <= (mp/2); v1++)
        for(int v2 = 0; v2 <= (mq/2); v2++)
        {
            mirp_G(Gy, fmp + mp, fmq + mq, mp, mq, v1, v2, gammap, gammaq, gammapq, working_prec);

            /* Gxy = Gx * Gy */
            arb_mul(Gxy, Gx, Gy, working_prec);

            for(int np = 0; np <= lmn1[2] + lmn2[2]; np++)
            for(int nq = 0; nq <= lmn3[2] + lmn4[2]; nq++)
            for(int w1 = 0; w1 <= (np/2); w1++)
            for(int w2 = 0; w2 <= (nq/2); w2++)
            {
                mirp_G(Gz, fnp + np, fnq + nq, np, nq, w1, w2, gammap, gammaq, gammapq, working_prec);

                /* Gxyz = Gx * Gy * Gz */
                arb_mul(Gxyz, Gxy, Gz, working_prec);

                for(int tx = 0; tx <= ((lp + lq - 2 * (u1 + u2)) / 2); tx++)
                for(int ty = 0; ty <= ((mp + mq - 2 * (v1 + v2)) / 2); ty++)
                for(int tz = 0; tz <= ((np + nq - 2 * (w1 + w2)) / 2); tz++)
                {
                    const int zeta = lp + lq + mp + mq + np + nq - 2*(u1 + u2 + v1 + v2 + w1 + w2) - tx - ty - tz;

                    const int xfac = lp + lq - 2*(u1 + u2 + tx);
                    const int yfac = mp + mq - 2*(v1 + v2 + ty);
                    const int zfac = np + nq - 2*(w1 + w2 + tz);


                    arb_set_si(tmp1, NEG1_POW(tx + ty + tz));
                    arb_mul(tmp1, tmp1, Gxyz, working_prec);

                    arb_mul(tmp1, tmp1, F + zeta, working_prec);

                    arb_pow_ui(tmp2, PQ+0, xfac, working_prec);
                    arb_mul(tmp1, tmp1, tmp2, working_prec);

                    arb_pow_ui(tmp2, PQ+1, yfac, working_prec);
                    arb_mul(tmp1, tmp1, tmp2, working_prec);

                    arb_pow_ui(tmp2, PQ+2, zfac, working_prec);
                    arb_mul(tmp1, tmp1, tmp2, working_prec);

                    arb_set_ui(tmp2, 4);
                    arb_pow_ui(tmp2, tmp2, u1 + u2 + tx + v1 + v2 + ty + w1 + w2 + tz, working_prec);

                    arb_pow_ui(tmp3, gammapq, tx + ty + tz, working_prec);
                    arb_mul(tmp2, tmp2, tmp3, working_prec);

                    arb_fac_ui(tmp3, xfac, working_prec);
                    arb_mul(tmp2, tmp2, tmp3, working_prec);

                    arb_fac_ui(tmp3, yfac, working_prec);
                    arb_mul(tmp2, tmp2, tmp3, working_prec);

                    arb_fac_ui(tmp3, zfac, working_prec);
                    arb_mul(tmp2, tmp2, tmp3, working_prec);

                    arb_fac_ui(tmp3, tx, working_prec);
                    arb_mul(tmp2, tmp2, tmp3, working_prec);

                    arb_fac_ui(tmp3, ty, working_prec);
                    arb_mul(tmp2, tmp2, tmp3, working_prec);

                    arb_fac_ui(tmp3, tz, working_prec);
                    arb_mul(tmp2, tmp2, tmp3, working_prec);

                    arb_div(tmp1, tmp1, tmp2, working_prec);

                    arb_add(integral, integral, tmp1, working_prec);
                }
            }
        }
    }


    /* Calculate the prefactor
     *
     * start with pfac = 2 * pi**2.5
     */
    arb_const_pi(tmp1, working_prec);
    arb_pow_ui(tmp1, tmp1, 5, working_prec);
    arb_sqrt(tmp1, tmp1, working_prec);
    arb_mul_ui(tmp1, tmp1, 2, working_prec);

    /*
     * Now multiply by K1 and K2
     * K1 = exp(-alpha1 * alpha2 * AB2 / gammap);
     * K2 = exp(-alpha3 * alpha4 * CD2 / gammaq);
     */
    arb_mul(tmp2, alpha1, alpha2, working_prec);
    arb_mul(tmp2, tmp2, AB2, working_prec);
    arb_div(tmp2, tmp2, gammap, working_prec);
    arb_mul_si(tmp2, tmp2, -1, working_prec);
    arb_exp(tmp2, tmp2, working_prec);
    arb_mul(tmp1, tmp1, tmp2, working_prec);

    arb_mul(tmp2, alpha3, alpha4, working_prec);
    arb_mul(tmp2, tmp2, CD2, working_prec);
    arb_div(tmp2, tmp2, gammaq, working_prec);
    arb_mul_si(tmp2, tmp2, -1, working_prec);
    arb_exp(tmp2, tmp2, working_prec);
    arb_mul(tmp1, tmp1, tmp2, working_prec);

    /*
     * divide by (gammap * gammaq * sqrt(gammap + gammaq))
     */
    arb_add(tmp2, gammap, gammaq, working_prec);
    arb_sqrt(tmp2, tmp2, working_prec);
    arb_mul(tmp2, tmp2, gammap, working_prec);
    arb_mul(tmp2, tmp2, gammaq, working_prec);
    arb_div(tmp1, tmp1, tmp2, working_prec);

    /* apply the prefactor */
    arb_mul(integral, integral, tmp1, working_prec);


    /* Check if this has no accuracy. This means the integral is zero */
    if(arb_rel_accuracy_bits(integral) < 0)
        arb_zero(integral);


    /* cleanup */
    _arb_vec_clear(F,   L+1);
    _arb_vec_clear(flp, lmn1[0]+lmn2[0]+1);
    _arb_vec_clear(fmp, lmn1[1]+lmn2[1]+1);
    _arb_vec_clear(fnp, lmn1[2]+lmn2[2]+1);
    _arb_vec_clear(flq, lmn3[0]+lmn4[0]+1);
    _arb_vec_clear(fmq, lmn3[1]+lmn4[1]+1);
    _arb_vec_clear(fnq, lmn3[2]+lmn4[2]+1);
    _arb_vec_clear(P,  3);
    _arb_vec_clear(PA, 3);
    _arb_vec_clear(PB, 3);
    _arb_vec_clear(Q,  3);
    _arb_vec_clear(QC, 3);
    _arb_vec_clear(QD, 3);
    _arb_vec_clear(PQ, 3);
    arb_clear(tmp1);
    arb_clear(tmp2);
    arb_clear(tmp3);
    arb_clear(gammap);
    arb_clear(gammaq);
    arb_clear(gammapq);
    arb_clear(AB2);
    arb_clear(CD2);
    arb_clear(PQ2);
    arb_clear(Gx);
    arb_clear(Gy);
    arb_clear(Gz);
    arb_clear(Gxy);
    arb_clear(Gxyz);
}

