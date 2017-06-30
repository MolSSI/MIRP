#include "mirp/kernels/boys.h"
#include "mirp/kernels/eri.h"
#include "mirp/arb_help.h"
#include "mirp/math.h"
#include "mirp/shell.h"
#include "mirp/gpt.h"

static void mirp_farr_interval(arb_t * f,
                               int lmn1, int lmn2,
                               arb_t xyz1, arb_t xyz2,
                               slong working_prec)
{
    int i, j, k;

    arb_t tmp1, tmp2;
    arb_init(tmp1);
    arb_init(tmp2);

    for (k = 0; k <= lmn1 + lmn2; k++)
    {
        arb_zero(f[k]);

        for (i = 0; i <= MIN(k,lmn1); i++)
        {
            j = k - i;
            if (j > lmn2)
                continue;

            mirp_binomial_coefficient_interval(tmp1, lmn1, i, working_prec);
            mirp_binomial_coefficient_interval(tmp2, lmn2, j, working_prec);
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
            arb_add(f[k], f[k], tmp1, working_prec);
        }
    }

    arb_clear(tmp1);
    arb_clear(tmp2);
}


static void mirp_G_interval(arb_t G, arb_t fp, arb_t fq,
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

void mirp_single_eri_interval(arb_t result,
                              const int * lmn1, const arb_t * A, const arb_t alpha1,
                              const int * lmn2, const arb_t * B, const arb_t alpha2,
                              const int * lmn3, const arb_t * C, const arb_t alpha3,
                              const int * lmn4, const arb_t * D, const arb_t alpha4,
                              slong working_prec)
{
    const int L_l = lmn1[0]+lmn2[0]+lmn3[0]+lmn4[0];
    const int L_m = lmn1[1]+lmn2[1]+lmn3[1]+lmn4[1];
    const int L_n = lmn1[2]+lmn2[2]+lmn3[2]+lmn4[2];
    const int L = L_l + L_m + L_n;

    arb_t F[L+1];
    arb_t flp[lmn1[0]+lmn2[0]+1];
    arb_t fmp[lmn1[1]+lmn2[1]+1];
    arb_t fnp[lmn1[2]+lmn2[2]+1];
    arb_t flq[lmn3[0]+lmn4[0]+1];
    arb_t fmq[lmn3[1]+lmn4[1]+1];
    arb_t fnq[lmn3[2]+lmn4[2]+1];
    mirp_init_arb_arr(F,   L+1    );
    mirp_init_arb_arr(flp, lmn1[0]+lmn2[0]+1);
    mirp_init_arb_arr(fmp, lmn1[1]+lmn2[1]+1);
    mirp_init_arb_arr(fnp, lmn1[2]+lmn2[2]+1);
    mirp_init_arb_arr(flq, lmn3[0]+lmn4[0]+1);
    mirp_init_arb_arr(fmq, lmn3[1]+lmn4[1]+1);
    mirp_init_arb_arr(fnq, lmn3[2]+lmn4[2]+1);

    /* Zero the result (we will be summing into it) */
    arb_zero(result);

    /* Temporary variables used in constructing expressions */
    arb_t tmp1, tmp2, tmp3;
    arb_init(tmp1);
    arb_init(tmp2);
    arb_init(tmp3);


    /*************************************************
     * Calculate all the various terms from the GPT
     *************************************************/
    arb_t gammap, P[3], PA[3], PB[3], AB2;
    arb_t gammaq, Q[3], QC[3], QD[3], CD2;
    arb_t gammapq, PQ[3], PQ2;
    mirp_init_arb_arr(P,  3);
    mirp_init_arb_arr(PA, 3);
    mirp_init_arb_arr(PB, 3);
    mirp_init_arb_arr(Q,  3);
    mirp_init_arb_arr(QC, 3);
    mirp_init_arb_arr(QD, 3);
    mirp_init_arb_arr(PQ, 3);
    arb_init(gammap);
    arb_init(gammaq);
    arb_init(gammapq);
    arb_init(AB2);
    arb_init(CD2);
    arb_init(PQ2);

    /* Gaussian Product Theorem */
    mirp_gpt_interval(alpha1, alpha2, A, B, gammap, P, PA, PB, AB2, working_prec);
    mirp_gpt_interval(alpha3, alpha4, C, D, gammaq, Q, QC, QD, CD2, working_prec);


    /*
     * gammapq = gammap * gammaq / (gammap + gammaq);
     * PQ[0] = P[0] - Q[0]
     * etc
     */
    arb_mul(tmp1,    gammap, gammaq, working_prec);
    arb_add(tmp2,    gammap, gammaq, working_prec);
    arb_div(gammapq, tmp1,   tmp2,   working_prec);

    arb_sub(PQ[0], P[0], Q[0], working_prec);
    arb_sub(PQ[1], P[1], Q[1], working_prec);
    arb_sub(PQ[2], P[2], Q[2], working_prec);

    /*
     * PQ2 = (P[0]-Q[0])*(P[0]-Q[0]) + (P[1]-Q[1])*(P[1]-Q[1]) + (P[2]-Q[2])*(P[2]-Q[2]);
     */
    arb_mul(PQ2, PQ[0], PQ[0], working_prec);
    arb_addmul(PQ2, PQ[1], PQ[1], working_prec);
    arb_addmul(PQ2, PQ[2], PQ[2], working_prec);


    mirp_farr_interval(flp, lmn1[0], lmn2[0], PA[0], PB[0], working_prec);
    mirp_farr_interval(fmp, lmn1[1], lmn2[1], PA[1], PB[1], working_prec);
    mirp_farr_interval(fnp, lmn1[2], lmn2[2], PA[2], PB[2], working_prec);
    mirp_farr_interval(flq, lmn3[0], lmn4[0], QC[0], QD[0], working_prec);
    mirp_farr_interval(fmq, lmn3[1], lmn4[1], QC[1], QD[1], working_prec);
    mirp_farr_interval(fnq, lmn3[2], lmn4[2], QC[2], QD[2], working_prec);


    /*
     *  Calculate the Boys function
     */
    arb_mul(tmp1, PQ2, gammapq, working_prec);
    mirp_boys_interval(F, L, tmp1, working_prec);


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
        mirp_G_interval(Gx, flp[lp], flq[lq], lp, lq, u1, u2, gammap, gammaq, gammapq, working_prec);

        for(int mp = 0; mp <= lmn1[1] + lmn2[1]; mp++)
        for(int mq = 0; mq <= lmn3[1] + lmn4[1]; mq++)
        for(int v1 = 0; v1 <= (mp/2); v1++)
        for(int v2 = 0; v2 <= (mq/2); v2++)
        {
            mirp_G_interval(Gy, fmp[mp], fmq[mq], mp, mq, v1, v2, gammap, gammaq, gammapq, working_prec);

            /* Gxy = Gx * Gy */
            arb_mul(Gxy, Gx, Gy, working_prec);

            for(int np = 0; np <= lmn1[2] + lmn2[2]; np++)
            for(int nq = 0; nq <= lmn3[2] + lmn4[2]; nq++)
            for(int w1 = 0; w1 <= (np/2); w1++)
            for(int w2 = 0; w2 <= (nq/2); w2++)
            {
                mirp_G_interval(Gz, fnp[np], fnq[nq], np, nq, w1, w2, gammap, gammaq, gammapq, working_prec);

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

                    arb_mul(tmp1, tmp1, F[zeta], working_prec);

                    arb_set_si(tmp2, xfac);
                    arb_pow(tmp2, PQ[0], tmp2, working_prec);
                    arb_mul(tmp1, tmp1, tmp2, working_prec);

                    arb_set_si(tmp2, yfac);
                    arb_pow(tmp2, PQ[1], tmp2, working_prec);
                    arb_mul(tmp1, tmp1, tmp2, working_prec);

                    arb_set_si(tmp2, zfac);
                    arb_pow(tmp2, PQ[2], tmp2, working_prec);
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

                    arb_add(result, result, tmp1, working_prec);
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
    arb_mul(result, result, tmp1, working_prec);


    /* cleanup */
    mirp_clear_arb_arr(F,   L+1);
    mirp_clear_arb_arr(flp, lmn1[0]+lmn2[0]+1);
    mirp_clear_arb_arr(fmp, lmn1[1]+lmn2[1]+1);
    mirp_clear_arb_arr(fnp, lmn1[2]+lmn2[2]+1);
    mirp_clear_arb_arr(flq, lmn3[0]+lmn4[0]+1);
    mirp_clear_arb_arr(fmq, lmn3[1]+lmn4[1]+1);
    mirp_clear_arb_arr(fnq, lmn3[2]+lmn4[2]+1);
    mirp_clear_arb_arr(P,  3);
    mirp_clear_arb_arr(PA, 3);
    mirp_clear_arb_arr(PB, 3);
    mirp_clear_arb_arr(Q,  3);
    mirp_clear_arb_arr(QC, 3);
    mirp_clear_arb_arr(QD, 3);
    mirp_clear_arb_arr(PQ, 3);
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


size_t mirp_prim_eri_interval(arb_t * result,
                              int am1, const arb_t * A, const arb_t alpha1,
                              int am2, const arb_t * B, const arb_t alpha2,
                              int am3, const arb_t * C, const arb_t alpha3,
                              int am4, const arb_t * D, const arb_t alpha4,
                              slong working_prec)
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
                    mirp_single_eri_interval(*(result + idx),
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


size_t mirp_eri_interval(arb_t * result,
                         int am1, const arb_t * A, int nprim1, int ngeneral1, const arb_t * alpha1, const arb_t * coeff1,
                         int am2, const arb_t * B, int nprim2, int ngeneral2, const arb_t * alpha2, const arb_t * coeff2,
                         int am3, const arb_t * C, int nprim3, int ngeneral3, const arb_t * alpha3, const arb_t * coeff3,
                         int am4, const arb_t * D, int nprim4, int ngeneral4, const arb_t * alpha4, const arb_t * coeff4,
                         slong working_prec)
{
    const size_t ncart1 = MIRP_NCART(am1);
    const size_t ncart2 = MIRP_NCART(am2);
    const size_t ncart3 = MIRP_NCART(am3);
    const size_t ncart4 = MIRP_NCART(am4);
    const size_t ncart1234 = ncart1*ncart2*ncart3*ncart4;
    const size_t ngeneral1234 = ngeneral1*ngeneral2*ngeneral3*ngeneral4;
    const size_t full_size = ncart1234*ngeneral1234;

    arb_t coeff;
    arb_init(coeff);

    arb_t * result_buffer = (arb_t *)malloc(full_size * sizeof(arb_t));
    mirp_init_arb_arr(result_buffer, full_size);

    for(size_t i = 0; i < full_size; i++)
        arb_zero(result[i]);

    for(int i = 0; i < nprim1; i++)
    for(int j = 0; j < nprim2; j++)
    for(int k = 0; k < nprim3; k++)
    for(int l = 0; l < nprim4; l++)
    {
        mirp_prim_eri_interval(result_buffer,
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
            arb_mul(coeff, coeff1[m*nprim1+i], coeff2[n*nprim2+j], working_prec);
            arb_mul(coeff, coeff,              coeff3[o*nprim3+k], working_prec);
            arb_mul(coeff, coeff,              coeff4[p*nprim4+l], working_prec);

            for(size_t q = 0; q < ncart1234; q++)
            {
                const size_t idx = ntotal*ncart1234+q;
                arb_addmul(result[idx], result_buffer[q], coeff, working_prec);
            }
            ntotal++;
        }
    }

    arb_clear(coeff);
    mirp_clear_arb_arr(result_buffer, full_size);
    free(result_buffer);    

    return full_size;
}
