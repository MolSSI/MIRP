#include "mirp/math.h"
#include "mirp/mpfr_help.h"
#include "mirp/boys.h"
#include "mirp/gpt.h"
#include <stdlib.h>

static void mirp_farr_mp(mpfr_t * f,
                           int lmn1, int lmn2,
                           mpfr_t xyz1, mpfr_t xyz2,
                           mpfr_prec_t working_prec)
{
    int i, j, k;

    mpfr_t tmp1, tmp2;
    mpfr_inits2(working_prec, tmp1, tmp2, (mpfr_ptr)0);

    for (k = 0; k <= lmn1 + lmn2; k++)
    {
        mpfr_set_zero(f[k], 0);

        for (i = 0; i <= MIN(k,lmn1); i++)
        {
            j = k - i;
            if (j > lmn2)
                continue;

            mirp_binomial_coefficient_mp(tmp1, lmn1, i);
            mirp_binomial_coefficient_mp(tmp2, lmn2, j);
            mpfr_mul(tmp1, tmp1, tmp2, MPFR_RNDN);

            if (lmn1 - i > 0)
            {
                mpfr_pow_ui(tmp2, xyz1, lmn1-i, MPFR_RNDN);
                mpfr_mul(tmp1, tmp1, tmp2, MPFR_RNDN);
            }
            if (lmn2 - j > 0)
            {
                mpfr_pow_ui(tmp2, xyz2, lmn2-j, MPFR_RNDN);
                mpfr_mul(tmp1, tmp1, tmp2, MPFR_RNDN);
            }
            mpfr_add(f[k], f[k], tmp1, MPFR_RNDN);
        }
    }

    mpfr_clear(tmp1);
    mpfr_clear(tmp2);
}


static void mirp_G_mp( mpfr_t G, mpfr_t fp, mpfr_t fq,
                        int np, int nq, int w1, int w2,
                        mpfr_t gammap, mpfr_t gammaq, mpfr_t gammapq,
                        mpfr_prec_t working_prec)
{
    mpfr_t tmp1, tmp2;
    mpfr_inits2(working_prec, tmp1, tmp2, (mpfr_ptr)0);

    mpfr_set_si(G, NEG1_POW(np), MPFR_RNDN);
    mpfr_mul(G, G, fp, MPFR_RNDN);
    mpfr_mul(G, G, fq, MPFR_RNDN);

    mpfr_fac_ui(tmp1, np, MPFR_RNDN);
    mpfr_mul(G, G, tmp1, MPFR_RNDN);
    mpfr_fac_ui(tmp1, nq, MPFR_RNDN);
    mpfr_mul(G, G, tmp1, MPFR_RNDN);

    mpfr_pow_si(tmp1, gammap, w1 - np, MPFR_RNDN);
    mpfr_mul(G, G, tmp1, MPFR_RNDN);
    mpfr_pow_si(tmp1, gammaq, w2 - nq, MPFR_RNDN);
    mpfr_mul(G, G, tmp1, MPFR_RNDN);

    mpfr_fac_ui(tmp1, np + nq - 2 * (w1 + w2), MPFR_RNDN);
    mpfr_mul(G, G, tmp1, MPFR_RNDN);

    mpfr_pow_ui(tmp1, gammapq, np + nq - 2 * (w1 + w2), MPFR_RNDN);
    mpfr_mul(G, G, tmp1, MPFR_RNDN);

    mpfr_fac_ui(tmp1, w1, MPFR_RNDN);
    mpfr_fac_ui(tmp2, w2, MPFR_RNDN);
    mpfr_mul(tmp1, tmp1, tmp2, MPFR_RNDN);
    mpfr_fac_ui(tmp2, np - 2 * w1, MPFR_RNDN);
    mpfr_mul(tmp1, tmp1, tmp2, MPFR_RNDN);
    mpfr_fac_ui(tmp2, nq - 2 * w2, MPFR_RNDN);
    mpfr_mul(tmp1, tmp1, tmp2, MPFR_RNDN);
    mpfr_div(G, G, tmp1, MPFR_RNDN);

    mpfr_clears(tmp1, tmp2, (mpfr_ptr)0);
}


void mirp_single_eri_mp(mpfr_t result,
                          int l1, int m1, int n1, mpfr_t alpha1, mpfr_t A[3],
                          int l2, int m2, int n2, mpfr_t alpha2, mpfr_t B[3],
                          int l3, int m3, int n3, mpfr_t alpha3, mpfr_t C[3],
                          int l4, int m4, int n4, mpfr_t alpha4, mpfr_t D[3],
                          mpfr_prec_t working_prec)
{
    const int L_l = l1+l2+l3+l4;
    const int L_m = m1+m2+m3+m4;
    const int L_n = n1+n2+n3+n4;
    const int L = L_l + L_m + L_n;

    mpfr_t F[L+1];
    mpfr_t flp[l1+l2+1];
    mpfr_t fmp[m1+m2+1];
    mpfr_t fnp[n1+n2+1];
    mpfr_t flq[l3+l4+1];
    mpfr_t fmq[m3+m4+1];
    mpfr_t fnq[n3+n4+1];
    mirp_init_mpfr_arr(F,   L+1,     working_prec);
    mirp_init_mpfr_arr(flp, l1+l2+1, working_prec);
    mirp_init_mpfr_arr(fmp, m1+m2+1, working_prec);
    mirp_init_mpfr_arr(fnp, n1+n2+1, working_prec);
    mirp_init_mpfr_arr(flq, l3+l4+1, working_prec);
    mirp_init_mpfr_arr(fmq, m3+m4+1, working_prec);
    mirp_init_mpfr_arr(fnq, n3+n4+1, working_prec);

    /* We need a temporary result variable in our working precision
     * We accumulate there, only (possibly) rounding at the very end
     */
    mpfr_t result_tmp;
    mpfr_init2(result_tmp, working_prec);
    mpfr_set_zero(result_tmp, 0);

    /* Temporary variables used in constructing expressions */
    mpfr_t tmp1, tmp2, tmp3;
    mpfr_inits2(working_prec, tmp1, tmp2, tmp3, (mpfr_ptr)0);


    /*************************************************
     * Calculate all the various terms from the GPT
     *************************************************/
    mpfr_t gammap, P[3], PA[3], PB[3], AB2;
    mpfr_t gammaq, Q[3], QC[3], QD[3], CD2;
    mpfr_t gammapq, PQ[3], PQ2;
    mirp_init_mpfr_arr(P,  3, working_prec);
    mirp_init_mpfr_arr(PA, 3, working_prec);
    mirp_init_mpfr_arr(PB, 3, working_prec);
    mirp_init_mpfr_arr(Q,  3, working_prec);
    mirp_init_mpfr_arr(QC, 3, working_prec);
    mirp_init_mpfr_arr(QD, 3, working_prec);
    mirp_init_mpfr_arr(PQ, 3, working_prec);
    mpfr_inits2(working_prec, gammap, gammaq, gammapq, AB2, CD2, PQ2, (mpfr_ptr)0);

    /* Gaussian Product Theorem */
    mirp_gpt_mp(alpha1, alpha2, A, B, gammap, P, PA, PB, AB2, working_prec);
    mirp_gpt_mp(alpha3, alpha4, C, D, gammaq, Q, QC, QD, CD2, working_prec);


    /*
     * gammapq = gammap * gammaq / (gammap + gammaq);
     * PQ[0] = P[0] - Q[0]
     * etc
     */
    mpfr_mul(tmp1,    gammap, gammaq, MPFR_RNDN);
    mpfr_add(tmp2,    gammap, gammaq, MPFR_RNDN);
    mpfr_div(gammapq, tmp1,   tmp2,   MPFR_RNDN);

    mpfr_sub(PQ[0], P[0], Q[0], MPFR_RNDN);
    mpfr_sub(PQ[1], P[1], Q[1], MPFR_RNDN);
    mpfr_sub(PQ[2], P[2], Q[2], MPFR_RNDN);

    /*
     * PQ2 = (P[0]-Q[0])*(P[0]-Q[0]) + (P[1]-Q[1])*(P[1]-Q[1]) + (P[2]-Q[2])*(P[2]-Q[2]);
     */
    mpfr_mul(PQ2, PQ[0], PQ[0], MPFR_RNDN);
    mpfr_fma(PQ2, PQ[1], PQ[1], PQ2, MPFR_RNDN);
    mpfr_fma(PQ2, PQ[2], PQ[2], PQ2, MPFR_RNDN);


    mirp_farr_mp(flp, l1, l2, PA[0], PB[0], working_prec);
    mirp_farr_mp(fmp, m1, m2, PA[1], PB[1], working_prec);
    mirp_farr_mp(fnp, n1, n2, PA[2], PB[2], working_prec);
    mirp_farr_mp(flq, l3, l4, QC[0], QD[0], working_prec);
    mirp_farr_mp(fmq, m3, m4, QC[1], QD[1], working_prec);
    mirp_farr_mp(fnq, n3, n4, QC[2], QD[2], working_prec);


    /*
     *  Calculate the Boys function
     */
    mpfr_mul(tmp1, PQ2, gammapq, MPFR_RNDN);
    mirp_boys_mp(F, L, tmp1, working_prec);


    /*
     * G values used within the loops
     */
    mpfr_t Gx, Gy, Gz, Gxy, Gxyz;
    mpfr_inits2(working_prec, Gx, Gy, Gz, Gxy, Gxyz, (mpfr_ptr)0);

    for(int lp = 0; lp <= l1 + l2; lp++)
    for(int lq = 0; lq <= l3 + l4; lq++)
    for(int u1 = 0; u1 <= (lp/2); u1++)
    for(int u2 = 0; u2 <= (lq/2); u2++)
    {
        mirp_G_mp(Gx, flp[lp], flq[lq], lp, lq, u1, u2, gammap, gammaq, gammapq, working_prec);

        for(int mp = 0; mp <= m1 + m2; mp++)
        for(int mq = 0; mq <= m3 + m4; mq++)
        for(int v1 = 0; v1 <= (mp/2); v1++)
        for(int v2 = 0; v2 <= (mq/2); v2++)
        {
            mirp_G_mp(Gy, fmp[mp], fmq[mq], mp, mq, v1, v2, gammap, gammaq, gammapq, working_prec);

            /* Gxy = Gx * Gy */
            mpfr_mul(Gxy, Gx, Gy, MPFR_RNDN);

            for(int np = 0; np <= n1 + n2; np++)
            for(int nq = 0; nq <= n3 + n4; nq++)
            for(int w1 = 0; w1 <= (np/2); w1++)
            for(int w2 = 0; w2 <= (nq/2); w2++)
            {
                mirp_G_mp(Gz, fnp[np], fnq[nq], np, nq, w1, w2, gammap, gammaq, gammapq, working_prec);

                /* Gxyz = Gx * Gy * Gz */
                mpfr_mul(Gxyz, Gxy, Gz, MPFR_RNDN);

                for(int tx = 0; tx <= ((lp + lq - 2 * (u1 + u2)) / 2); tx++)
                for(int ty = 0; ty <= ((mp + mq - 2 * (v1 + v2)) / 2); ty++)
                for(int tz = 0; tz <= ((np + nq - 2 * (w1 + w2)) / 2); tz++)
                {
                    const int zeta = lp + lq + mp + mq + np + nq - 2*(u1 + u2 + v1 + v2 + w1 + w2) - tx - ty - tz;

                    const int xfac = lp + lq - 2*(u1 + u2 + tx);
                    const int yfac = mp + mq - 2*(v1 + v2 + ty);
                    const int zfac = np + nq - 2*(w1 + w2 + tz);


                    mpfr_set_si(tmp1, NEG1_POW(tx + ty + tz), MPFR_RNDN);
                    mpfr_mul(tmp1, tmp1, Gxyz, MPFR_RNDN);

                    mpfr_mul(tmp1, tmp1, F[zeta], MPFR_RNDN);

                    mpfr_pow_si(tmp2, PQ[0], xfac, MPFR_RNDN);
                    mpfr_mul(tmp1, tmp1, tmp2, MPFR_RNDN);

                    mpfr_pow_si(tmp2, PQ[1], yfac, MPFR_RNDN);
                    mpfr_mul(tmp1, tmp1, tmp2, MPFR_RNDN);

                    mpfr_pow_si(tmp2, PQ[2], zfac, MPFR_RNDN);
                    mpfr_mul(tmp1, tmp1, tmp2, MPFR_RNDN);

                    mpfr_set_ui(tmp2, 4, MPFR_RNDN);
                    mpfr_pow_si(tmp2, tmp2, u1 + u2 + tx + v1 + v2 + ty + w1 + w2 + tz, MPFR_RNDN);

                    mpfr_pow_si(tmp3, gammapq, tx + ty + tz, MPFR_RNDN);
                    mpfr_mul(tmp2, tmp2, tmp3, MPFR_RNDN);

                    mpfr_fac_ui(tmp3, xfac, MPFR_RNDN);
                    mpfr_mul(tmp2, tmp2, tmp3, MPFR_RNDN);

                    mpfr_fac_ui(tmp3, yfac, MPFR_RNDN);
                    mpfr_mul(tmp2, tmp2, tmp3, MPFR_RNDN);

                    mpfr_fac_ui(tmp3, zfac, MPFR_RNDN);
                    mpfr_mul(tmp2, tmp2, tmp3, MPFR_RNDN);

                    mpfr_fac_ui(tmp3, tx, MPFR_RNDN);
                    mpfr_mul(tmp2, tmp2, tmp3, MPFR_RNDN);

                    mpfr_fac_ui(tmp3, ty, MPFR_RNDN);
                    mpfr_mul(tmp2, tmp2, tmp3, MPFR_RNDN);

                    mpfr_fac_ui(tmp3, tz, MPFR_RNDN);
                    mpfr_mul(tmp2, tmp2, tmp3, MPFR_RNDN);

                    mpfr_div(tmp1, tmp1, tmp2, MPFR_RNDN);

                    mpfr_add(result_tmp, result_tmp, tmp1, MPFR_RNDN);
                }
            }
        }
    }


    /* Calculate the prefactor
     *
     * start with pfac = 2 * pi**2.5
     */
    mpfr_const_pi(tmp1, MPFR_RNDN);
    mpfr_pow_ui(tmp1, tmp1, 5, MPFR_RNDN);
    mpfr_sqrt(tmp1, tmp1, MPFR_RNDN);
    mpfr_mul_ui(tmp1, tmp1, 2, MPFR_RNDN);

    /*
     * Now multiply by K1 and K2
     * K1 = exp(-alpha1 * alpha2 * AB2 / gammap);
     * K2 = exp(-alpha3 * alpha4 * CD2 / gammaq);
     */
    mpfr_mul(tmp2, alpha1, alpha2, MPFR_RNDN);
    mpfr_mul(tmp2, tmp2, AB2, MPFR_RNDN);
    mpfr_div(tmp2, tmp2, gammap, MPFR_RNDN);
    mpfr_setsign(tmp2, tmp2, 1, MPFR_RNDN);
    mpfr_exp(tmp2, tmp2, MPFR_RNDN);
    mpfr_mul(tmp1, tmp1, tmp2, MPFR_RNDN);

    mpfr_mul(tmp2, alpha3, alpha4, MPFR_RNDN);
    mpfr_mul(tmp2, tmp2, CD2, MPFR_RNDN);
    mpfr_div(tmp2, tmp2, gammaq, MPFR_RNDN);
    mpfr_setsign(tmp2, tmp2, 1, MPFR_RNDN);
    mpfr_exp(tmp2, tmp2, MPFR_RNDN);
    mpfr_mul(tmp1, tmp1, tmp2, MPFR_RNDN);

    /*
     * divide by (gammap * gammaq * sqrt(gammap + gammaq))
     */
    mpfr_add(tmp2, gammap, gammaq, MPFR_RNDN);
    mpfr_sqrt(tmp2, tmp2, MPFR_RNDN);
    mpfr_mul(tmp2, tmp2, gammap, MPFR_RNDN);
    mpfr_mul(tmp2, tmp2, gammaq, MPFR_RNDN);
    mpfr_div(tmp1, tmp1, tmp2, MPFR_RNDN);

    /* apply the prefactor */
    mpfr_mul(result_tmp, result_tmp, tmp1, MPFR_RNDN);

    /* Store in the actual result area passed to this function (possibly rounding) */
    mpfr_set(result, result_tmp, MPFR_RNDN);

    /* cleanup */
    mirp_clear_mpfr_arr(F,   L+1);
    mirp_clear_mpfr_arr(flp, l1+l2+1);
    mirp_clear_mpfr_arr(fmp, m1+m2+1);
    mirp_clear_mpfr_arr(fnp, n1+n2+1);
    mirp_clear_mpfr_arr(flq, l3+l4+1);
    mirp_clear_mpfr_arr(fmq, m3+m4+1);
    mirp_clear_mpfr_arr(fnq, n3+n4+1);
    mirp_clear_mpfr_arr(P,  3);
    mirp_clear_mpfr_arr(PA, 3);
    mirp_clear_mpfr_arr(PB, 3);
    mirp_clear_mpfr_arr(Q,  3);
    mirp_clear_mpfr_arr(QC, 3);
    mirp_clear_mpfr_arr(QD, 3);
    mirp_clear_mpfr_arr(PQ, 3);
    mpfr_clears(tmp1, tmp2, tmp3, (mpfr_ptr)0);
    mpfr_clears(gammap, gammaq, gammapq, AB2, CD2, PQ2, (mpfr_ptr)0);
    mpfr_clears(Gx, Gy, Gz, Gxy, Gxyz, (mpfr_ptr)0);
    mpfr_clear(result_tmp);
}


void mirp_single_eri_mp_str(char ** result,
                              int l1, int m1, int n1, const char * alpha1, const char * A[3],
                              int l2, int m2, int n2, const char * alpha2, const char * B[3],
                              int l3, int m3, int n3, const char * alpha3, const char * C[3],
                              int l4, int m4, int n4, const char * alpha4, const char * D[3],
                              long working_prec)
{
    mpfr_t alpha1_mp, alpha2_mp, alpha3_mp, alpha4_mp;
    mpfr_t A_mp[3], B_mp[3], C_mp[3], D_mp[3];
    mpfr_t result_mp;
    mpfr_init2(result_mp, (mpfr_prec_t)working_prec);

    mpfr_init_set_str(alpha1_mp, alpha1, 10, MPFR_RNDN);
    mpfr_init_set_str(alpha2_mp, alpha2, 10, MPFR_RNDN);
    mpfr_init_set_str(alpha3_mp, alpha3, 10, MPFR_RNDN);
    mpfr_init_set_str(alpha4_mp, alpha4, 10, MPFR_RNDN);
    mpfr_init_set_str(A_mp[0], A[0], 10, MPFR_RNDN);
    mpfr_init_set_str(A_mp[1], A[1], 10, MPFR_RNDN);
    mpfr_init_set_str(A_mp[2], A[2], 10, MPFR_RNDN);
    mpfr_init_set_str(B_mp[0], B[0], 10, MPFR_RNDN);
    mpfr_init_set_str(B_mp[1], B[1], 10, MPFR_RNDN);
    mpfr_init_set_str(B_mp[2], B[2], 10, MPFR_RNDN);
    mpfr_init_set_str(C_mp[0], C[0], 10, MPFR_RNDN);
    mpfr_init_set_str(C_mp[1], C[1], 10, MPFR_RNDN);
    mpfr_init_set_str(C_mp[2], C[2], 10, MPFR_RNDN);
    mpfr_init_set_str(D_mp[0], D[0], 10, MPFR_RNDN);
    mpfr_init_set_str(D_mp[1], D[1], 10, MPFR_RNDN);
    mpfr_init_set_str(D_mp[2], D[2], 10, MPFR_RNDN);

    mirp_single_eri_mp(result_mp,
                         l1, m1, n1, alpha1_mp, A_mp,
                         l2, m2, n2, alpha2_mp, B_mp,
                         l3, m3, n3, alpha3_mp, C_mp,
                         l4, m4, n4, alpha4_mp, D_mp,
                         (mpfr_prec_t)working_prec);

   
    mpfr_asprintf(result, "%Re", result_mp);

    mirp_clear_mpfr_arr(A_mp, 3);
    mirp_clear_mpfr_arr(B_mp, 3);
    mirp_clear_mpfr_arr(C_mp, 3);
    mirp_clear_mpfr_arr(D_mp, 3);
    mpfr_clears(alpha1_mp, alpha2_mp, alpha3_mp, alpha4_mp, (mpfr_ptr)0);
    mpfr_clear(result_mp);
}
