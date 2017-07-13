/*! \file
 *
 * \brief Calculation of electron repulsion integrals (double precision)
 */

#include "mirp/kernels/boys.h"
#include "mirp/kernels/eri.h"
#include "mirp/loops.h"
#include "mirp/math.h"
#include "mirp/gpt.h"
#include <math.h>


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

void mirp_single_eri_double(double * output,
                            const int * lmn1, const double * A, double alpha1,
                            const int * lmn2, const double * B, double alpha2,
                            const int * lmn3, const double * C, double alpha3,
                            const int * lmn4, const double * D, double alpha4)
{
    const int L_l = lmn1[0]+lmn2[0]+lmn3[0]+lmn4[0];
    const int L_m = lmn1[1]+lmn2[1]+lmn3[1]+lmn4[1];
    const int L_n = lmn1[2]+lmn2[2]+lmn3[2]+lmn4[2];
    const int L = L_l + L_m + L_n;

    *output = 0.0;

    double F[L+1];
    double flp[lmn1[0]+lmn2[0]+1];
    double fmp[lmn1[1]+lmn2[1]+1];
    double fnp[lmn1[2]+lmn2[2]+1];
    double flq[lmn3[0]+lmn4[0]+1];
    double fmq[lmn3[1]+lmn4[1]+1];
    double fnq[lmn3[2]+lmn4[2]+1];

    double gammap, P[3], PA[3], PB[3], AB2;
    double gammaq, Q[3], QC[3], QD[3], CD2;

    mirp_gpt_double(alpha1, alpha2, A, B, &gammap, P, PA, PB, &AB2);
    mirp_gpt_double(alpha3, alpha4, C, D, &gammaq, Q, QC, QD, &CD2);

    const double gammapq = gammap * gammaq / (gammap + gammaq);
    const double PQ[3] = { P[0] - Q[0],
                           P[1] - Q[1],
                           P[2] - Q[2] };

    const double PQ2 = PQ[0]*PQ[0] + PQ[1]*PQ[1] + PQ[2]*PQ[2];

    compute_farr(flp, lmn1[0], lmn2[0], PA[0], PB[0]);
    compute_farr(fmp, lmn1[1], lmn2[1], PA[1], PB[1]);
    compute_farr(fnp, lmn1[2], lmn2[2], PA[2], PB[2]);
    compute_farr(flq, lmn3[0], lmn4[0], QC[0], QD[0]);
    compute_farr(fmq, lmn3[1], lmn4[1], QC[1], QD[1]);
    compute_farr(fnq, lmn3[2], lmn4[2], QC[2], QD[2]);

    mirp_boys_double(F, L, PQ2 * gammapq);

    for(int lp = 0; lp <= lmn1[0] + lmn2[0]; lp++)
    for(int lq = 0; lq <= lmn3[0] + lmn4[0]; lq++)
    for(int u1 = 0; u1 <= (lp/2); u1++)
    for(int u2 = 0; u2 <= (lq/2); u2++)
    {
        double Gx = NEG1_POW(lp) * flp[lp] * flq[lq] * mirp_factorial(lp) * mirp_factorial(lq)
                  * pow(gammap, u1 - lp) * pow(gammaq, u2 - lq) * mirp_factorial(lp + lq - 2 * (u1 + u2))
                  * pow(gammapq, lp + lq - 2 * (u1 + u2));
        Gx /= (mirp_factorial(u1) * mirp_factorial(u2) * mirp_factorial(lp - 2 * u1) * mirp_factorial(lq - 2 * u2));

        for(int mp = 0; mp <= lmn1[1] + lmn2[1]; mp++)
        for(int mq = 0; mq <= lmn3[1] + lmn4[1]; mq++)
        for(int v1 = 0; v1 <= (mp/2); v1++)
        for(int v2 = 0; v2 <= (mq/2); v2++)
        {
            double Gy = NEG1_POW(mp) * fmp[mp] * fmq[mq] * mirp_factorial(mp) * mirp_factorial(mq)
                      * pow(gammap, v1 - mp) * pow(gammaq, v2 - mq) * mirp_factorial(mp + mq - 2 * (v1 + v2))
                      * pow(gammapq, mp + mq - 2 * (v1 + v2));
            Gy /= (mirp_factorial(v1) * mirp_factorial(v2) * mirp_factorial(mp - 2 * v1) * mirp_factorial(mq - 2 * v2));

            for(int np = 0; np <= lmn1[2] + lmn2[2]; np++)
            for(int nq = 0; nq <= lmn3[2] + lmn4[2]; nq++)
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

                    *output += tmp;
                }
            }
        }
    }

    double pfac = 2 * pow(MIRP_PI, 2.5)
                * exp(-alpha1 * alpha2 * AB2 / gammap)
                * exp(-alpha3 * alpha4 * CD2 / gammaq);
    pfac /= (gammap * gammaq * sqrt(gammap + gammaq));

    *output *= pfac;
}

void mirp_prim_eri_double(double * output,
                          int am1, const double * A, double alpha1,
                          int am2, const double * B, double alpha2,
                          int am3, const double * C, double alpha3,
                          int am4, const double * D, double alpha4)
{
    mirp_cartloop4_double(output,
                          am1, A, alpha1,
                          am2, B, alpha2,
                          am3, C, alpha3,
                          am4, D, alpha4,
                          mirp_single_eri_double); // callback
}


void mirp_eri_double(double * output,
                     int am1, const double * A, int nprim1, int ngeneral1, const double * alpha1, const double * coeff1,
                     int am2, const double * B, int nprim2, int ngeneral2, const double * alpha2, const double * coeff2,
                     int am3, const double * C, int nprim3, int ngeneral3, const double * alpha3, const double * coeff3,
                     int am4, const double * D, int nprim4, int ngeneral4, const double * alpha4, const double * coeff4)
{
    mirp_loop4_double(output,
                      am1, A, nprim1, ngeneral1, alpha1, coeff1,
                      am2, B, nprim2, ngeneral2, alpha2, coeff2,
                      am3, C, nprim3, ngeneral3, alpha3, coeff3,
                      am4, D, nprim4, ngeneral4, alpha4, coeff4,
                      mirp_prim_eri_double);
}

