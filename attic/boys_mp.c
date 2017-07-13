/*! \file
 *
 * \brief Calculation of the boys function using arbitrary precision
 */

#include "mirp/kernels/boys.h"
#include "mirp/mpfr_help.h"
#include <assert.h>

void mirp_boys_mp(mpfr_t *F, int m, const mpfr_t t, mpfr_prec_t working_prec)
{
    assert(m >= 0);
    assert(working_prec > 0);

    int i;

    mpfr_t t2, et, sum, term, lastterm, test;
    mpfr_t tmp1, tmp2, tmp3;
    mpfr_inits2(working_prec, t2, et, sum, term, lastterm, test, (mpfr_ptr)0);
    mpfr_inits2(working_prec, tmp1, tmp2, tmp3, (mpfr_ptr)0);

    mpfr_t Ftmp[m+1]; /* in the working precision */
    mirp_init_mpfr_arr(Ftmp, m+1, working_prec);

    /* t2 = 2*t */
    mpfr_mul_si(t2, t, 2, MPFR_RNDN);

    /* et = exp(-t)
       Note: t is always positive, so we can use mpfr_neg
     */
    mpfr_neg(et, t, MPFR_RNDN);
    mpfr_exp(et, et, MPFR_RNDN);

    int do_short = 0;


    /* The short-range formula converges much better for
     * t < (m + 3/2)
     * So skip the long range if that happens
     * Note that this is a conservative bound, and could probably be increased
     */
    mpfr_set_si(test, 2*m+3, MPFR_RNDN);
    mpfr_div_si(test, test, 2, MPFR_RNDN);
         
    if(mpfr_less_p(t, test))
        do_short = 1;
    
    if(!do_short)
    {
        /* Attempt the long-range approximation */ 
        mpfr_const_pi(tmp1, MPFR_RNDN);
        mpfr_div(tmp1, tmp1, t, MPFR_RNDN);
        mpfr_sqrt(tmp1, tmp1, MPFR_RNDN);
        mpfr_div_si(tmp1, tmp1, 2, MPFR_RNDN);
        
        for(i = 1; i <= m; i++)
        {
            mpfr_set_si(tmp2, 2*i-1, MPFR_RNDN);
            mpfr_div(tmp2, tmp2, t2, MPFR_RNDN);
            mpfr_mul(tmp1, tmp1, tmp2, MPFR_RNDN);
        }
        mpfr_set(Ftmp[m], tmp1, MPFR_RNDN);

        /* Determine the error associated with the long-range approximation */
        mpfr_set_zero(sum, 0);
        mpfr_set_si(term, 1, MPFR_RNDN);

        i = 0;
        do {
            i++;
            mpfr_abs(lastterm, term, MPFR_RNDN);
            mpfr_mul_si(term, term, 2*m - 2*i + 1, MPFR_RNDN);
            mpfr_div(term, term, t2, MPFR_RNDN);
            mpfr_abs(test, term, MPFR_RNDN);

            if(mpfr_greater_p(test, lastterm))
                break;

            mpfr_set(test, sum, MPFR_RNDN);
            mpfr_add(sum, sum, term, MPFR_RNDN);
            mpfr_sub(test, test, sum, MPFR_RNDN);

        } while(!mpfr_zero_p(test));

        mpfr_mul(sum, sum, et, MPFR_RNDN);
        mpfr_div(sum, sum, t2, MPFR_RNDN);

        /*
         * Determine if this error is satisfactory
         * If not, mark that we have to do the short-range version
         */
        mpfr_add(test, Ftmp[m], sum, MPFR_RNDN);
        if(!mpfr_equal_p(test, Ftmp[m]))
            do_short = 1;
    }

    if(do_short)
    {
        /* The long-range approximation isn't good enough */
        mpfr_set_si(sum, 1, MPFR_RNDN);
        mpfr_set_si(term, 1, MPFR_RNDN);

        i = 0;
        do {
            i++;
            mpfr_mul(term, term, t2, MPFR_RNDN);
            mpfr_div_si(term, term, 2*m + 2*i + 1, MPFR_RNDN);

            /* store the old term, then update and calculate the difference */
            mpfr_set(test, sum, MPFR_RNDN);
            mpfr_add(sum, sum, term, MPFR_RNDN);
            mpfr_sub(test, test, sum, MPFR_RNDN);
        } while (!mpfr_zero_p(test));

        mpfr_mul(sum, sum, et, MPFR_RNDN);
        mpfr_div_si(Ftmp[m], sum, 2*m+1, MPFR_RNDN);
    }

    /* Now do downwards recursion */
    mpfr_set(F[m], Ftmp[m], MPFR_RNDN);
    for(i = m - 1; i >= 0; i--)
    {
        /* F[i] = (t2 * F[i + 1] + et) / (2 * i + 1) */
        mpfr_fma(Ftmp[i], t2, Ftmp[i + 1], et, MPFR_RNDN);
        mpfr_div_si(Ftmp[i], Ftmp[i], 2 * i + 1, MPFR_RNDN);

        /* move the result into the original array (with
         * its original precision */
        mpfr_set(F[i], Ftmp[i], MPFR_RNDN);
    }

    mpfr_clears(t2, et, sum, term, lastterm, test, (mpfr_ptr)0);
    mpfr_clears(tmp1, tmp2, tmp3, (mpfr_ptr)0);
    mirp_clear_mpfr_arr(Ftmp, m+1);
}

