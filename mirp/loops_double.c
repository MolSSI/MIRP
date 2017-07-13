/*! \file
 *
 * \brief Helpers for looping over cartesian functions and primitives
 */

#include "mirp/loops.h"
#include "mirp/shell.h"
#include <string.h> /* for memset */

void mirp_cartloop4_double(double * output,
                           int am1, const double * A, double alpha1,
                           int am2, const double * B, double alpha2,
                           int am3, const double * C, double alpha3,
                           int am4, const double * D, double alpha4,
                           cb_single4_double cb)
{
    const long ncart1 = MIRP_NCART(am1);
    const long ncart2 = MIRP_NCART(am2);
    const long ncart3 = MIRP_NCART(am3);
    const long ncart4 = MIRP_NCART(am4);

    long idx = 0;
    int lmn1[3] = {am1, 0, 0};
    for(long i = 0; i < ncart1; i++)
    {
        int lmn2[3] = {am2, 0, 0};
        for(long j = 0; j < ncart2; j++)
        {
            int lmn3[3] = {am3, 0, 0};
            for(long k = 0; k < ncart3; k++)
            {
                int lmn4[3] = {am4, 0, 0};
                for(long l = 0; l < ncart4; l++)
                {
                    cb(output + idx,
                       lmn1, A, alpha1,
                       lmn2, B, alpha2,
                       lmn3, C, alpha3,
                       lmn4, D, alpha4);

                    idx++;

                    mirp_iterate_gaussian(lmn4);
                }

                mirp_iterate_gaussian(lmn3);
            }

            mirp_iterate_gaussian(lmn2);
        }

        mirp_iterate_gaussian(lmn1);
    }
}


void mirp_loop4_double(double * output,
                       int am1, const double * A, int nprim1, int ngeneral1, const double * alpha1, const double * coeff1,
                       int am2, const double * B, int nprim2, int ngeneral2, const double * alpha2, const double * coeff2,
                       int am3, const double * C, int nprim3, int ngeneral3, const double * alpha3, const double * coeff3,
                       int am4, const double * D, int nprim4, int ngeneral4, const double * alpha4, const double * coeff4,
                       cb_prim4_double cb)
{
    const long ncart1 = MIRP_NCART(am1);
    const long ncart2 = MIRP_NCART(am2);
    const long ncart3 = MIRP_NCART(am3);
    const long ncart4 = MIRP_NCART(am4);
    const long ncart1234 = ncart1*ncart2*ncart3*ncart4;
    const long ngeneral1234 = ngeneral1*ngeneral2*ngeneral3*ngeneral4;
    const long full_size = ncart1234*ngeneral1234;

    double * output_buffer = malloc(full_size * sizeof(double));
    double * coeff1_norm = malloc(nprim1 * ngeneral1 * sizeof(double));
    double * coeff2_norm = malloc(nprim2 * ngeneral2 * sizeof(double));
    double * coeff3_norm = malloc(nprim3 * ngeneral3 * sizeof(double));
    double * coeff4_norm = malloc(nprim4 * ngeneral4 * sizeof(double));

    mirp_normalize_shell_double(am1, nprim1, ngeneral1, alpha1, coeff1, coeff1_norm);
    mirp_normalize_shell_double(am2, nprim2, ngeneral2, alpha2, coeff2, coeff2_norm);
    mirp_normalize_shell_double(am3, nprim3, ngeneral3, alpha3, coeff3, coeff3_norm);
    mirp_normalize_shell_double(am4, nprim4, ngeneral4, alpha4, coeff4, coeff4_norm);

    memset(output, 0, full_size * sizeof(double));
    for(int i = 0; i < nprim1; i++)
    for(int j = 0; j < nprim2; j++)
    for(int k = 0; k < nprim3; k++)
    for(int l = 0; l < nprim4; l++)
    {
        cb(output_buffer,
           am1, A, alpha1[i],
           am2, B, alpha2[j],
           am3, C, alpha3[k],
           am4, D, alpha4[l]);

        long ntotal = 0;
        for(int m = 0; m < ngeneral1; m++)
        for(int n = 0; n < ngeneral2; n++)
        for(int o = 0; o < ngeneral3; o++)
        for(int p = 0; p < ngeneral4; p++)
        {
            const double coeff = coeff1_norm[m*nprim1+i]
                               * coeff2_norm[n*nprim2+j]
                               * coeff3_norm[o*nprim3+k]
                               * coeff4_norm[p*nprim4+l];

            for(long q = 0; q < ncart1234; q++)
                output[ntotal*ncart1234+q] += output_buffer[q] * coeff;
            ntotal++;
        }
    }

    free(output_buffer);    
    free(coeff1_norm);
    free(coeff2_norm);
    free(coeff3_norm);
    free(coeff4_norm);
}

