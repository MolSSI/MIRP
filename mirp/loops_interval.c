/*! \file
 *
 * \brief Helpers for looping over cartesian functions and primitives
 */

#include "mirp/loops.h"
#include "mirp/shell.h"

/*! \brief Compute all cartesian components of a single primitive integral
 *         (interval arithmetic)
 *
 * \copydetails mirp_cartloop4_double
 * \param [in] working_prec The working precision (binary digits/bits) to use in the calculation
 */
static void mirp_cartloop4_interval(arb_ptr output,
                                    int am1, arb_srcptr A, const arb_t alpha1,
                                    int am2, arb_srcptr B, const arb_t alpha2,
                                    int am3, arb_srcptr C, const arb_t alpha3,
                                    int am4, arb_srcptr D, const arb_t alpha4,
                                    slong working_prec, cb_single4_interval cb)
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


void mirp_loop4_interval(arb_ptr output,
                        int am1, arb_srcptr A, int nprim1, int ngeneral1, arb_srcptr alpha1, arb_srcptr coeff1,
                        int am2, arb_srcptr B, int nprim2, int ngeneral2, arb_srcptr alpha2, arb_srcptr coeff2,
                        int am3, arb_srcptr C, int nprim3, int ngeneral3, arb_srcptr alpha3, arb_srcptr coeff3,
                        int am4, arb_srcptr D, int nprim4, int ngeneral4, arb_srcptr alpha4, arb_srcptr coeff4,
                        slong working_prec, cb_single4_interval cb)
{
    const long ncart1 = MIRP_NCART(am1);
    const long ncart2 = MIRP_NCART(am2);
    const long ncart3 = MIRP_NCART(am3);
    const long ncart4 = MIRP_NCART(am4);
    const long ncart1234 = ncart1*ncart2*ncart3*ncart4;
    const long ngeneral1234 = ngeneral1*ngeneral2*ngeneral3*ngeneral4;
    const long full_size = ncart1234*ngeneral1234;

    arb_ptr output_buffer = _arb_vec_init(full_size);
    arb_ptr coeff1_norm = _arb_vec_init(nprim1 * ngeneral1);
    arb_ptr coeff2_norm = _arb_vec_init(nprim2 * ngeneral2);
    arb_ptr coeff3_norm = _arb_vec_init(nprim3 * ngeneral3);
    arb_ptr coeff4_norm = _arb_vec_init(nprim4 * ngeneral4);

    mirp_normalize_shell_interval(am1, nprim1, ngeneral1, alpha1, coeff1, coeff1_norm, working_prec);
    mirp_normalize_shell_interval(am2, nprim2, ngeneral2, alpha2, coeff2, coeff2_norm, working_prec);
    mirp_normalize_shell_interval(am3, nprim3, ngeneral3, alpha3, coeff3, coeff3_norm, working_prec);
    mirp_normalize_shell_interval(am4, nprim4, ngeneral4, alpha4, coeff4, coeff4_norm, working_prec);

    _arb_vec_zero(output, full_size);

    /* A temporary variable (used to build up the coefficient) */
    arb_t coeff;
    arb_init(coeff);

    for(int i = 0; i < nprim1; i++)
    for(int j = 0; j < nprim2; j++)
    for(int k = 0; k < nprim3; k++)
    for(int l = 0; l < nprim4; l++)
    {
        mirp_cartloop4_interval(output_buffer,
                                am1, A, alpha1 + i,
                                am2, B, alpha2 + j,
                                am3, C, alpha3 + k,
                                am4, D, alpha4 + l,
                                working_prec, cb);

        long ntotal = 0;
        for(int m = 0; m < ngeneral1; m++)
        for(int n = 0; n < ngeneral2; n++)
        for(int o = 0; o < ngeneral3; o++)
        for(int p = 0; p < ngeneral4; p++)
        {
            arb_mul(coeff, coeff1_norm+(m*nprim1+i), coeff2_norm+(n*nprim2+j), working_prec);
            arb_mul(coeff, coeff,                    coeff3_norm+(o*nprim3+k), working_prec);
            arb_mul(coeff, coeff,                    coeff4_norm+(p*nprim4+l), working_prec);

            for(long q = 0; q < ncart1234; q++)
            {
                const long idx = ntotal*ncart1234 + q;
                arb_addmul(output+idx, output_buffer+q, coeff, working_prec);
            }
            ntotal++;
        }
    }

    arb_clear(coeff);
    _arb_vec_clear(output_buffer, full_size);
    _arb_vec_clear(coeff1_norm, nprim1 * ngeneral1);
    _arb_vec_clear(coeff2_norm, nprim2 * ngeneral2);
    _arb_vec_clear(coeff3_norm, nprim3 * ngeneral3);
    _arb_vec_clear(coeff4_norm, nprim4 * ngeneral4);
}

