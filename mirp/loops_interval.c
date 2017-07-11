/*! \file
 *
 * \brief Helpers for looping over cartesian functions and primitives
 */

#include "mirp/loops.h"
#include "mirp/shell.h"
#include "mirp/math.h"
#include "mirp/arb_help.h"

void mirp_prim_cartloop4_interval(arb_t * result,
                                  int am1, const arb_t * A, const arb_t alpha1,
                                  int am2, const arb_t * B, const arb_t alpha2,
                                  int am3, const arb_t * C, const arb_t alpha3,
                                  int am4, const arb_t * D, const arb_t alpha4,
                                  mpfr_prec_t working_prec, cb_4gaussians_interval cb)
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
                    cb(result[idx],
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
}

void mirp_cartloop4_interval(arb_t * result,
                            int am1, const arb_t * A, int nprim1, int ngeneral1, const arb_t * alpha1, const arb_t * coeff1,
                            int am2, const arb_t * B, int nprim2, int ngeneral2, const arb_t * alpha2, const arb_t * coeff2,
                            int am3, const arb_t * C, int nprim3, int ngeneral3, const arb_t * alpha3, const arb_t * coeff3,
                            int am4, const arb_t * D, int nprim4, int ngeneral4, const arb_t * alpha4, const arb_t * coeff4,
                            mpfr_prec_t working_prec, cb_4shells_interval cb)
{
    const long ncart1 = MIRP_NCART(am1);
    const long ncart2 = MIRP_NCART(am2);
    const long ncart3 = MIRP_NCART(am3);
    const long ncart4 = MIRP_NCART(am4);
    const long ncart1234 = ncart1*ncart2*ncart3*ncart4;
    const long ngeneral1234 = ngeneral1*ngeneral2*ngeneral3*ngeneral4;
    const long full_size = ncart1234*ngeneral1234;

    arb_t coeff;
    arb_init(coeff);

    arb_t * result_buffer = (arb_t *)malloc(full_size * sizeof(arb_t));
    arb_t * coeff1_norm = malloc(nprim1 * ngeneral1 * sizeof(arb_t));
    arb_t * coeff2_norm = malloc(nprim2 * ngeneral2 * sizeof(arb_t));
    arb_t * coeff3_norm = malloc(nprim3 * ngeneral3 * sizeof(arb_t));
    arb_t * coeff4_norm = malloc(nprim4 * ngeneral4 * sizeof(arb_t));

    mirp_init_arb_arr(result_buffer, full_size);
    mirp_init_arb_arr(coeff1_norm, nprim1 * ngeneral1);
    mirp_init_arb_arr(coeff2_norm, nprim2 * ngeneral2);
    mirp_init_arb_arr(coeff3_norm, nprim3 * ngeneral3);
    mirp_init_arb_arr(coeff4_norm, nprim4 * ngeneral4);

    mirp_normalize_shell_interval(am1, nprim1, ngeneral1, alpha1, coeff1, coeff1_norm, working_prec);
    mirp_normalize_shell_interval(am2, nprim2, ngeneral2, alpha2, coeff2, coeff2_norm, working_prec);
    mirp_normalize_shell_interval(am3, nprim3, ngeneral3, alpha3, coeff3, coeff3_norm, working_prec);
    mirp_normalize_shell_interval(am4, nprim4, ngeneral4, alpha4, coeff4, coeff4_norm, working_prec);


    for(long i = 0; i < full_size; i++)
        arb_zero(result[i]);

    for(int i = 0; i < nprim1; i++)
    for(int j = 0; j < nprim2; j++)
    for(int k = 0; k < nprim3; k++)
    for(int l = 0; l < nprim4; l++)
    {
        cb(result_buffer,
           am1, A, alpha1[i],
           am2, B, alpha2[j],
           am3, C, alpha3[k],
           am4, D, alpha4[l],
           working_prec);

        long ntotal = 0;
        for(int m = 0; m < ngeneral1; m++)
        for(int n = 0; n < ngeneral2; n++)
        for(int o = 0; o < ngeneral3; o++)
        for(int p = 0; p < ngeneral4; p++)
        {
            arb_mul(coeff, coeff1_norm[m*nprim1+i], coeff2_norm[n*nprim2+j], working_prec);
            arb_mul(coeff, coeff,                   coeff3_norm[o*nprim3+k], working_prec);
            arb_mul(coeff, coeff,                   coeff4_norm[p*nprim4+l], working_prec);

            for(long q = 0; q < ncart1234; q++)
            {
                const long idx = ntotal*ncart1234+q;
                arb_addmul(result[idx], result_buffer[q], coeff, working_prec);
            }
            ntotal++;
        }
    }

    arb_clear(coeff);
    mirp_clear_arb_arr(result_buffer, full_size);
    mirp_clear_arb_arr(coeff1_norm, nprim1 * ngeneral1);
    mirp_clear_arb_arr(coeff2_norm, nprim2 * ngeneral2);
    mirp_clear_arb_arr(coeff3_norm, nprim3 * ngeneral3);
    mirp_clear_arb_arr(coeff4_norm, nprim4 * ngeneral4);
    free(result_buffer);    
    free(coeff1_norm);
    free(coeff2_norm);
    free(coeff3_norm);
    free(coeff4_norm);
}

