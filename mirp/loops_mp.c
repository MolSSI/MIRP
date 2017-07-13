/*! \file
 *
 * \brief Helpers for looping over cartesian functions and primitives
 */

#include "mirp/loops.h"
#include "mirp/shell.h"
#include "mirp/math.h"
#include "mirp/mpfr_help.h"

void mirp_cartloop4_mp(mpfr_t * result,
                       int am1, const mpfr_t * A, const mpfr_t alpha1,
                       int am2, const mpfr_t * B, const mpfr_t alpha2,
                       int am3, const mpfr_t * C, const mpfr_t alpha3,
                       int am4, const mpfr_t * D, const mpfr_t alpha4,
                       mpfr_prec_t working_prec, cb_single4_mp cb)
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

void mirp_loop4_mp(mpfr_t * result,
                   int am1, const mpfr_t * A, int nprim1, int ngeneral1, const mpfr_t * alpha1, const mpfr_t * coeff1,
                   int am2, const mpfr_t * B, int nprim2, int ngeneral2, const mpfr_t * alpha2, const mpfr_t * coeff2,
                   int am3, const mpfr_t * C, int nprim3, int ngeneral3, const mpfr_t * alpha3, const mpfr_t * coeff3,
                   int am4, const mpfr_t * D, int nprim4, int ngeneral4, const mpfr_t * alpha4, const mpfr_t * coeff4,
                   mpfr_prec_t working_prec, cb_prim4_mp cb)
{

    const long ncart1 = MIRP_NCART(am1);
    const long ncart2 = MIRP_NCART(am2);
    const long ncart3 = MIRP_NCART(am3);
    const long ncart4 = MIRP_NCART(am4);
    const long ncart1234 = ncart1*ncart2*ncart3*ncart4;
    const long ngeneral1234 = ngeneral1*ngeneral2*ngeneral3*ngeneral4;
    const long full_size = ncart1234*ngeneral1234;

    mpfr_t coeff;
    mpfr_init2(coeff, working_prec);

    mpfr_t * result_tmp = (mpfr_t *)malloc(full_size * sizeof(mpfr_t));
    mpfr_t * result_buffer = (mpfr_t *)malloc(full_size * sizeof(mpfr_t));
    mirp_init_mpfr_arr(result_tmp, full_size, working_prec);
    mirp_init_mpfr_arr(result_buffer, full_size, working_prec);

    for(long i = 0; i < full_size; i++)
        mpfr_set_zero(result_tmp[i], 0);

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
            mpfr_mul(coeff, coeff1[m*nprim1+i], coeff2[n*nprim2+j], MPFR_RNDN);
            mpfr_mul(coeff, coeff,              coeff3[o*nprim3+k], MPFR_RNDN);
            mpfr_mul(coeff, coeff,              coeff4[p*nprim4+l], MPFR_RNDN);

            for(long q = 0; q < ncart1234; q++)
            {
                const long idx = ntotal*ncart1234+q;
                mpfr_fma(result_tmp[idx], result_buffer[q], coeff, result_tmp[idx], MPFR_RNDN);
            }
            ntotal++;
        }
    }

    for(long i = 0; i < full_size; i++)
        mpfr_set(result[i], result_tmp[i], MPFR_RNDN);

    mpfr_clear(coeff);
    mirp_clear_mpfr_arr(result_tmp, full_size);
    mirp_clear_mpfr_arr(result_buffer, full_size);
    free(result_buffer);    
    free(result_tmp);    
}

