#include <string.h> /* for memset */
#include <math.h>
#include "mirp/math.h"
#include "mirp/shell.h"
#include "mirp/boys.h"
#include "mirp/gpt.h"
#include "mirp/eri.h"


static void compute_farr(double * f, int lmn1, int lmn2, double xyz1, double xyz2)
{
    int i, j, k;
    double tmp;

    for (k = 0; k <= lmn1 + lmn2; k++)
    {
        f[k] = 0.0;
        for (i = 0; i <= MIN(k,lmn1); i++)
        {
            j = k - i;
            if (j > lmn2)
                continue;
            tmp = mirp_binomial_coefficient(lmn1,i) * mirp_binomial_coefficient(lmn2,j);
            if (lmn1 - i > 0)
                tmp *= pow(xyz1, lmn1 - i);
            if (lmn2 - j > 0)
                tmp *= pow(xyz2, lmn2 - j);
            f[k] += tmp;
        }
    }
}

void mirp_single_eri_double(double * result,
                            int l1, int m1, int n1, double alpha1, const double A[3],
                            int l2, int m2, int n2, double alpha2, const double B[3],
                            int l3, int m3, int n3, double alpha3, const double C[3],
                            int l4, int m4, int n4, double alpha4, const double D[3])
{
    const int L_l = l1+l2+l3+l4;
    const int L_m = m1+m2+m3+m4;
    const int L_n = n1+n2+n3+n4;
    const int L = L_l + L_m + L_n;

    *result = 0.0;

    double F[L+1];
    double flp[l1+l2+1];
    double fmp[m1+m2+1];
    double fnp[n1+n2+1];
    double flq[l3+l4+1];
    double fmq[m3+m4+1];
    double fnq[n3+n4+1];

    double gammap, P[3], PA[3], PB[3], AB2;
    double gammaq, Q[3], QC[3], QD[3], CD2;

    mirp_gpt(alpha1, alpha2, A, B, &gammap, P, PA, PB, &AB2);
    mirp_gpt(alpha3, alpha4, C, D, &gammaq, Q, QC, QD, &CD2);

    const double gammapq = gammap * gammaq / (gammap + gammaq);
    const double PQ[3] = { P[0] - Q[0],
                           P[1] - Q[1],
                           P[2] - Q[2] };

    const double PQ2 = PQ[0]*PQ[0] + PQ[1]*PQ[1] + PQ[2]*PQ[2];

    compute_farr(flp, l1, l2, PA[0], PB[0]);
    compute_farr(fmp, m1, m2, PA[1], PB[1]);
    compute_farr(fnp, n1, n2, PA[2], PB[2]);
    compute_farr(flq, l3, l4, QC[0], QD[0]);
    compute_farr(fmq, m3, m4, QC[1], QD[1]);
    compute_farr(fnq, n3, n4, QC[2], QD[2]);

    mirp_boys_double(F, L, PQ2 * gammapq);

    for(int lp = 0; lp <= l1 + l2; lp++)
    for(int lq = 0; lq <= l3 + l4; lq++)
    for(int u1 = 0; u1 <= (lp/2); u1++)
    for(int u2 = 0; u2 <= (lq/2); u2++)
    {
        double Gx = NEG1_POW(lp) * flp[lp] * flq[lq] * mirp_factorial(lp) * mirp_factorial(lq)
                  * pow(gammap, u1 - lp) * pow(gammaq, u2 - lq) * mirp_factorial(lp + lq - 2 * (u1 + u2))
                  * pow(gammapq, lp + lq - 2 * (u1 + u2));
        Gx /= (mirp_factorial(u1) * mirp_factorial(u2) * mirp_factorial(lp - 2 * u1) * mirp_factorial(lq - 2 * u2));

        for(int mp = 0; mp <= m1 + m2; mp++)
        for(int mq = 0; mq <= m3 + m4; mq++)
        for(int v1 = 0; v1 <= (mp/2); v1++)
        for(int v2 = 0; v2 <= (mq/2); v2++)
        {
            double Gy = NEG1_POW(mp) * fmp[mp] * fmq[mq] * mirp_factorial(mp) * mirp_factorial(mq)
                      * pow(gammap, v1 - mp) * pow(gammaq, v2 - mq) * mirp_factorial(mp + mq - 2 * (v1 + v2))
                      * pow(gammapq, mp + mq - 2 * (v1 + v2));
            Gy /= (mirp_factorial(v1) * mirp_factorial(v2) * mirp_factorial(mp - 2 * v1) * mirp_factorial(mq - 2 * v2));

            for(int np = 0; np <= n1 + n2; np++)
            for(int nq = 0; nq <= n3 + n4; nq++)
            for(int w1 = 0; w1 <= (np/2); w1++)
            for(int w2 = 0; w2 <= (nq/2); w2++)
            {
                double Gz = NEG1_POW(np) * fnp[np] * fnq[nq] * mirp_factorial(np) * mirp_factorial(nq)
                          * pow(gammap, w1 - np) * pow(gammaq, w2 - nq) * mirp_factorial(np + nq - 2 * (w1 + w2))
                          * pow(gammapq, np + nq - 2 * (w1 + w2));
                Gz /= (mirp_factorial(w1) * mirp_factorial(w2) * mirp_factorial(np - 2 * w1) * mirp_factorial(nq - 2 * w2));

                double Gxyz = Gx*Gy*Gz;

                for(int tx = 0; tx <= ((lp + lq - 2 * (u1 + u2)) / 2); tx++)
                for(int ty = 0; ty <= ((mp + mq - 2 * (v1 + v2)) / 2); ty++)
                for(int tz = 0; tz <= ((np + nq - 2 * (w1 + w2)) / 2); tz++)
                {
                    const int zeta = lp + lq + mp + mq + np + nq - 2*(u1 + u2 + v1 + v2 + w1 + w2) - tx - ty - tz;

                    const int xfac = lp + lq - 2*(u1 + u2 + tx);
                    const int yfac = mp + mq - 2*(v1 + v2 + ty);
                    const int zfac = np + nq - 2*(w1 + w2 + tz);

                    double tmp = NEG1_POW(tx + ty + tz) * Gxyz * F[zeta]
                               * pow(PQ[0], xfac) * pow(PQ[1], yfac) * pow(PQ[2], zfac);

                    tmp /= pow(4, u1 + u2 + tx + v1 + v2 + ty + w1 + w2 + tz)
                         * pow(gammapq, tx + ty + tz)
                         * mirp_factorial(xfac) * mirp_factorial(tx)
                         * mirp_factorial(yfac) * mirp_factorial(ty)
                         * mirp_factorial(zfac) * mirp_factorial(tz);

                    *result += tmp;
                }
            }
        }
    }

    double pfac = 2 * pow(MIRP_PI, 2.5)
                * exp(-alpha1 * alpha2 * AB2 / gammap)
                * exp(-alpha3 * alpha4 * CD2 / gammaq);
    pfac /= (gammap * gammaq * sqrt(gammap + gammaq));

    *result *= pfac;
}

size_t mirp_prim_eri_double(double * result,
                            int am1, double alpha1, const double A[3],
                            int am2, double alpha2, const double B[3],
                            int am3, double alpha3, const double C[3],
                            int am4, double alpha4, const double D[3])
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
                    mirp_single_eri_double(result + idx,
                                           lmn1[0], lmn1[1], lmn1[2], alpha1, A,
                                           lmn2[0], lmn2[1], lmn2[2], alpha2, B,
                                           lmn3[0], lmn3[1], lmn3[2], alpha3, C,
                                           lmn4[0], lmn4[1], lmn4[2], alpha4, D);

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


size_t mirp_eri_double(double * result,
                       int am1, const double A[3], int nprim1, int ngeneral1, const double * alpha1, const double * coeff1,
                       int am2, const double B[3], int nprim2, int ngeneral2, const double * alpha2, const double * coeff2,
                       int am3, const double C[3], int nprim3, int ngeneral3, const double * alpha3, const double * coeff3,
                       int am4, const double D[4], int nprim4, int ngeneral4, const double * alpha4, const double * coeff4)
{

    const size_t ncart1 = MIRP_NCART(am1);
    const size_t ncart2 = MIRP_NCART(am2);
    const size_t ncart3 = MIRP_NCART(am3);
    const size_t ncart4 = MIRP_NCART(am4);
    const size_t ncart1234 = ncart1*ncart2*ncart3*ncart4;
    const size_t ngeneral1234 = ngeneral1*ngeneral2*ngeneral3*ngeneral4;
    const size_t full_size = ncart1234*ngeneral1234;

    double * result_buffer = (double *)malloc(full_size * sizeof(double));
    memset(result, 0, full_size * sizeof(double));

    for(int i = 0; i < nprim1; i++)
    for(int j = 0; j < nprim2; j++)
    for(int k = 0; k < nprim3; k++)
    for(int l = 0; l < nprim4; l++)
    {
        mirp_prim_eri_double(result_buffer,
                             am1, alpha1[i], A,
                             am2, alpha2[j], B,
                             am3, alpha3[k], C,
                             am4, alpha4[l], D);

        size_t ntotal = 0;
        for(int m = 0; m < ngeneral1; m++)
        for(int n = 0; n < ngeneral2; n++)
        for(int o = 0; o < ngeneral3; o++)
        for(int p = 0; p < ngeneral4; p++)
        {
            const double coeff = coeff1[m*nprim1+i]
                               * coeff2[n*nprim2+j]
                               * coeff3[o*nprim3+k]
                               * coeff4[p*nprim4+l];

            for(size_t q = 0; q < ncart1234; q++)
                result[ntotal*ncart1234+q] += result_buffer[q] * coeff;
            ntotal++;
        }

    }

    free(result_buffer);    

    return full_size;
}

