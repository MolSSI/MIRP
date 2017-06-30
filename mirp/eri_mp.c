#include "mirp/math.h"
#include "mirp/shell.h"
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
    mpfr_init2(tmp1, working_prec);
    mpfr_init2(tmp2, working_prec);

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


static void mirp_G_mp(mpfr_t G, mpfr_t fp, mpfr_t fq,
                      int np, int nq, int w1, int w2,
                      mpfr_t gammap, mpfr_t gammaq, mpfr_t gammapq,
                      mpfr_prec_t working_prec)
{
    mpfr_t tmp1, tmp2;
    mpfr_init2(tmp1, working_prec);
    mpfr_init2(tmp2, working_prec);

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
                        const int * lmn1, const mpfr_t * A, const mpfr_t alpha1,
                        const int * lmn2, const mpfr_t * B, const mpfr_t alpha2,
                        const int * lmn3, const mpfr_t * C, const mpfr_t alpha3,
                        const int * lmn4, const mpfr_t * D, const mpfr_t alpha4,
                        mpfr_prec_t working_prec)
{
    const int L_l = lmn1[0]+lmn2[0]+lmn3[0]+lmn4[0];
    const int L_m = lmn1[1]+lmn2[1]+lmn3[1]+lmn4[1];
    const int L_n = lmn1[2]+lmn2[2]+lmn3[2]+lmn4[2];
    const int L = L_l + L_m + L_n;

    mpfr_t F[L+1];
    mpfr_t flp[lmn1[0]+lmn2[0]+1];
    mpfr_t fmp[lmn1[1]+lmn2[1]+1];
    mpfr_t fnp[lmn1[2]+lmn2[2]+1];
    mpfr_t flq[lmn3[0]+lmn4[0]+1];
    mpfr_t fmq[lmn3[1]+lmn4[1]+1];
    mpfr_t fnq[lmn3[2]+lmn4[2]+1];
    mirp_init_mpfr_arr(F,   L+1,     working_prec);
    mirp_init_mpfr_arr(flp, lmn1[0]+lmn2[0]+1, working_prec);
    mirp_init_mpfr_arr(fmp, lmn1[1]+lmn2[1]+1, working_prec);
    mirp_init_mpfr_arr(fnp, lmn1[2]+lmn2[2]+1, working_prec);
    mirp_init_mpfr_arr(flq, lmn3[0]+lmn4[0]+1, working_prec);
    mirp_init_mpfr_arr(fmq, lmn3[1]+lmn4[1]+1, working_prec);
    mirp_init_mpfr_arr(fnq, lmn3[2]+lmn4[2]+1, working_prec);

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


    mirp_farr_mp(flp, lmn1[0], lmn2[0], PA[0], PB[0], working_prec);
    mirp_farr_mp(fmp, lmn1[1], lmn2[1], PA[1], PB[1], working_prec);
    mirp_farr_mp(fnp, lmn1[2], lmn2[2], PA[2], PB[2], working_prec);
    mirp_farr_mp(flq, lmn3[0], lmn4[0], QC[0], QD[0], working_prec);
    mirp_farr_mp(fmq, lmn3[1], lmn4[1], QC[1], QD[1], working_prec);
    mirp_farr_mp(fnq, lmn3[2], lmn4[2], QC[2], QD[2], working_prec);


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

    for(int lp = 0; lp <= lmn1[0] + lmn2[0]; lp++)
    for(int lq = 0; lq <= lmn3[0] + lmn4[0]; lq++)
    for(int u1 = 0; u1 <= (lp/2); u1++)
    for(int u2 = 0; u2 <= (lq/2); u2++)
    {
        mirp_G_mp(Gx, flp[lp], flq[lq], lp, lq, u1, u2, gammap, gammaq, gammapq, working_prec);

        for(int mp = 0; mp <= lmn1[1] + lmn2[1]; mp++)
        for(int mq = 0; mq <= lmn3[1] + lmn4[1]; mq++)
        for(int v1 = 0; v1 <= (mp/2); v1++)
        for(int v2 = 0; v2 <= (mq/2); v2++)
        {
            mirp_G_mp(Gy, fmp[mp], fmq[mq], mp, mq, v1, v2, gammap, gammaq, gammapq, working_prec);

            /* Gxy = Gx * Gy */
            mpfr_mul(Gxy, Gx, Gy, MPFR_RNDN);

            for(int np = 0; np <= lmn1[2] + lmn2[2]; np++)
            for(int nq = 0; nq <= lmn3[2] + lmn4[2]; nq++)
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

                    mpfr_pow_ui(tmp3, gammapq, tx + ty + tz, MPFR_RNDN);
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
    mirp_clear_mpfr_arr(flp, lmn1[0]+lmn2[0]+1);
    mirp_clear_mpfr_arr(fmp, lmn1[1]+lmn2[1]+1);
    mirp_clear_mpfr_arr(fnp, lmn1[2]+lmn2[2]+1);
    mirp_clear_mpfr_arr(flq, lmn3[0]+lmn4[0]+1);
    mirp_clear_mpfr_arr(fmq, lmn3[1]+lmn4[1]+1);
    mirp_clear_mpfr_arr(fnq, lmn3[2]+lmn4[2]+1);
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


size_t mirp_prim_eri_mp(mpfr_t * result,
                        int am1, const mpfr_t * A, const mpfr_t alpha1,
                        int am2, const mpfr_t * B, const mpfr_t alpha2,
                        int am3, const mpfr_t * C, const mpfr_t alpha3,
                        int am4, const mpfr_t * D, const mpfr_t alpha4,
                        mpfr_prec_t working_prec)
{

    const size_t ncart1 = MIRP_NCART(am1);
    const size_t ncart2 = MIRP_NCART(am2);
    const size_t ncart3 = MIRP_NCART(am3);
    const size_t ncart4 = MIRP_NCART(am4);
    const size_t ncart1234 = ncart1*ncart2*ncart3*ncart4;

    size_t idx = 0;
    int lmn1[3] = {am1, 0, 0};
    for(size_t i = 0; i < ncart1; i++)
    {
        int lmn2[3] = {am2, 0, 0};
        for(size_t j = 0; j < ncart2; j++)
        {
            int lmn3[3] = {am3, 0, 0};
            for(size_t k = 0; k < ncart3; k++)
            {
                int lmn4[3] = {am4, 0, 0};
                for(size_t l = 0; l < ncart4; l++)
                {
                    mirp_single_eri_mp(*(result + idx),
                                       lmn1, A, alpha1,
                                       lmn2, B, alpha2,
                                       lmn3, C, alpha3,
                                       lmn4, D, alpha4,
                                       working_prec);

                    idx++;

                    mirp_iterate_gaussian(lmn4);
                }

                mirp_iterate_gaussian(lmn3);
            }

            mirp_iterate_gaussian(lmn2);
        }

        mirp_iterate_gaussian(lmn1);
    }

    return ncart1234;
}

size_t mirp_eri_mp(mpfr_t * result,
                   int am1, const mpfr_t * A, int nprim1, int ngeneral1, const mpfr_t * alpha1, const mpfr_t * coeff1,
                   int am2, const mpfr_t * B, int nprim2, int ngeneral2, const mpfr_t * alpha2, const mpfr_t * coeff2,
                   int am3, const mpfr_t * C, int nprim3, int ngeneral3, const mpfr_t * alpha3, const mpfr_t * coeff3,
                   int am4, const mpfr_t * D, int nprim4, int ngeneral4, const mpfr_t * alpha4, const mpfr_t * coeff4,
                   mpfr_prec_t working_prec)
{
    const size_t ncart1 = MIRP_NCART(am1);
    const size_t ncart2 = MIRP_NCART(am2);
    const size_t ncart3 = MIRP_NCART(am3);
    const size_t ncart4 = MIRP_NCART(am4);
    const size_t ncart1234 = ncart1*ncart2*ncart3*ncart4;
    const size_t ngeneral1234 = ngeneral1*ngeneral2*ngeneral3*ngeneral4;
    const size_t full_size = ncart1234*ngeneral1234;

    mpfr_t coeff;
    mpfr_init2(coeff, working_prec);

    mpfr_t * result_tmp = (mpfr_t *)malloc(full_size * sizeof(mpfr_t));
    mpfr_t * result_buffer = (mpfr_t *)malloc(full_size * sizeof(mpfr_t));
    mirp_init_mpfr_arr(result_tmp, full_size, working_prec);
    mirp_init_mpfr_arr(result_buffer, full_size, working_prec);

    for(size_t i = 0; i < full_size; i++)
        mpfr_set_zero(result_tmp[i], 0);

    for(int i = 0; i < nprim1; i++)
    for(int j = 0; j < nprim2; j++)
    for(int k = 0; k < nprim3; k++)
    for(int l = 0; l < nprim4; l++)
    {
        mirp_prim_eri_mp(result_buffer,
                         am1, A, alpha1[i],
                         am2, B, alpha2[j],
                         am3, C, alpha3[k],
                         am4, D, alpha4[l],
                         working_prec);


        size_t ntotal = 0;
        for(int m = 0; m < ngeneral1; m++)
        for(int n = 0; n < ngeneral2; n++)
        for(int o = 0; o < ngeneral3; o++)
        for(int p = 0; p < ngeneral4; p++)
        {
            mpfr_mul(coeff, coeff1[m*nprim1+i], coeff2[n*nprim2+j], MPFR_RNDN);
            mpfr_mul(coeff, coeff,              coeff3[o*nprim3+k], MPFR_RNDN);
            mpfr_mul(coeff, coeff,              coeff4[p*nprim4+l], MPFR_RNDN);

            for(size_t q = 0; q < ncart1234; q++)
            {
                const size_t idx = ntotal*ncart1234+q;
                mpfr_fma(result_tmp[idx], result_buffer[q], coeff, result_tmp[idx], MPFR_RNDN);
            }
            ntotal++;
        }
    }

    for(size_t i = 0; i < full_size; i++)
        mpfr_set(result[i], result_tmp[i], MPFR_RNDN);

    mpfr_clear(coeff);
    mirp_clear_mpfr_arr(result_tmp, full_size);
    mirp_clear_mpfr_arr(result_buffer, full_size);
    free(result_buffer);    
    free(result_tmp);    

    return full_size;
}
