/*! \file
 *
 * \brief Calculation of the boys function in double precision
 *        and arbitrary precision
 */

#include "mirp/boys.h"
#include "mirp/math.h"
#include "mirp/arb_help.h"
#include <assert.h>

void mirp_boys_interval(arb_t *F, int m, const arb_t t, slong working_prec)
{
    assert(m >= 0);
    assert(working_prec > 0);

    int i;

    arb_t t2, et, sum, term, lastterm, test;
    arb_t tmp1, tmp2, tmp3;
    arb_init(t2);
    arb_init(et);
    arb_init(sum);
    arb_init(term);
    arb_init(lastterm);
    arb_init(test);
    arb_init(tmp1);
    arb_init(tmp2);
    arb_init(tmp3);

    /* t2 = 2*t */
    arb_mul_ui(t2, t, 2, working_prec);

    /* et = exp(-x)
       Note: x is always positive, so we can use arb_neg
     */
    arb_neg(et, t);
    arb_exp(et, et, working_prec);

    int do_short = 0;

    /* The short-range formula converges much better for
     * t < (m + 3/2)
     * So skip the long range if that happens
     * Note that this is a conservative bound, and could probably be increased
     */
    arb_set_si(test, 2*m+3);
    arb_div_si(test, test, 2, working_prec);

    if(arb_lt(t, test))
        do_short = 1;

    if(!do_short)
    {

        /* Value of the long-range approximation */
        mirp_double_factorial_interval(tmp1, 2*m-1, working_prec);
        arb_set_ui(tmp2, 2);
        arb_pow_ui(tmp2, tmp2, m+1, working_prec);
        arb_div(tmp1, tmp1, tmp2, working_prec);

        arb_const_pi(tmp2, working_prec);
        arb_pow_ui(tmp3, t, 2*m+1, working_prec);
        arb_div(tmp2, tmp2, tmp3, working_prec);
        arb_sqrt(tmp2, tmp2, working_prec);
        arb_mul(F[m], tmp1, tmp2, working_prec);

        /* Determine the error associated with the long-range approximation */
        arb_zero(sum);
        arb_set_ui(term, 1);

        i = 0;
        do {
            i++;
            arb_abs(lastterm, term);
            arb_mul_si(term, term, 2*m - 2*i + 1, working_prec);
            arb_div(term, term, t2, working_prec);
            arb_abs(test, term);

            if(arb_gt(test, lastterm))
                break;

            arb_set(test, sum);
            arb_add(sum, sum, term, working_prec);
        } while(!arb_contains(sum, test)); /* Is the old term contained completely
                                              within the new term */

        arb_mul(sum, sum, et, working_prec);
        arb_div(sum, sum, t2, working_prec);

        /*
         * Determine if this error is satisfactory
         * If not, mark that we have to do the short-range version
         */
        if(!arb_overlaps(F[m], sum))
            do_short = 1;
    }

    if(do_short)
    {
        arb_set_ui(sum, 1);
        arb_set_ui(term, 1);
    
        i = 0;
        do
        {
            i++;
            arb_mul(term, term, t2, working_prec);
            arb_div_si(term, term, 2*m + 2*i + 1, working_prec);

            /* store the old term, then update and calculate the difference */
            arb_set(test, sum);
            arb_add(sum, sum, term, working_prec);
        } while(!arb_contains(sum, test)); /* Is the old term contained completely
                                              within the new term */

        arb_mul(F[m], sum, et, working_prec);
        arb_div_si(F[m], F[m], 2*m+1, working_prec);
    }
        
    /* Now do downwards recursion */
    for(i = m - 1; i >= 0; i--)
    {
        /* F[m] = (t2 * F[m + 1] + et) / (2 * m + 1) */
        arb_mul(F[i], t2, F[i + 1], working_prec);
        arb_add(F[i], F[i], et, working_prec);
        arb_div_ui(F[i], F[i], 2 * i + 1, working_prec);
    }

    arb_init(t2);
    arb_init(et);
    arb_init(sum);
    arb_init(term);
    arb_init(lastterm);
    arb_init(test);
    arb_init(tmp1);
    arb_init(tmp2);
    arb_init(tmp3);
}

