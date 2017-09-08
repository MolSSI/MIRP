/*! \file
 *
 * \brief Calculation of the boys function using interval arithmetic
 */

#include "mirp/kernels/boys.h"
#include "mirp/math.h"
#include <assert.h>

void mirp_boys(arb_ptr F, int m, const arb_t t, slong working_prec)
{
    assert(m >= 0);
    assert(!(arb_is_negative(t)));
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
        /* Attempt the long-range approximation */
        arb_const_pi(tmp1, working_prec);
        arb_div(tmp1, tmp1, t, working_prec);
        arb_sqrt(tmp1, tmp1, working_prec);
        arb_div_si(tmp1, tmp1, 2, working_prec);

        for(i = 1; i <= m; i++)
        {
            arb_set_si(tmp2, 2*i-1);
            arb_div(tmp2, tmp2, t2, working_prec);
            arb_mul(tmp1, tmp1, tmp2, working_prec);
        }

        arb_set(F + m, tmp1);

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
        } while(!arb_contains(sum, test)); /* Is the old term contained
                                              completely within the new term */

        /* Note about the above test: If the old term is contained within
           the new term, all we have done is added error. Ie, we've
           added such a small number that the midpoint hasn't changed (much),
           but the error has increased */

        //printf("Done with long-range test in %d cycles\n", i);

        arb_mul(sum, sum, et, working_prec);
        arb_div(sum, sum, t2, working_prec);

        /*
         * Determine if this error is satisfactory
         * If not, mark that we have to do the short-range version
         */
        arb_sub(test, F+m, sum, working_prec);
        if(!arb_overlaps(F+m, test))
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
        } while(!arb_contains(sum, test)); /* Is the old term contained
                                              completely within the new term */

        arb_mul(F+m, sum, et, working_prec);
        arb_div_si(F+m, F+m, 2*m+1, working_prec);
        //printf("Done with short-range approximation in %d cycles\n", i);
    }

    /* Now do downwards recursion */
    for(i = m - 1; i >= 0; i--)
    {
        /* F+m = (t2 * F[m + 1] + et) / (2 * m + 1) */
        arb_mul(F+i, t2, F + (i + 1), working_prec);
        arb_add(F+i, F+i, et, working_prec);
        arb_div_ui(F+i, F+i, 2 * i + 1, working_prec);
    }

    arb_clear(t2);
    arb_clear(et);
    arb_clear(sum);
    arb_clear(term);
    arb_clear(lastterm);
    arb_clear(test);
    arb_clear(tmp1);
    arb_clear(tmp2);
    arb_clear(tmp3);
}


void mirp_boys_str(arb_ptr F, int m, const char * t, slong working_prec)
{
    arb_t t_mp;
    arb_init(t_mp);

    arb_set_str(t_mp, t, working_prec);
    mirp_boys(F, m, t_mp, working_prec);

    arb_clear(t_mp);
}


void mirp_boys_exact(double *F, int m, double t)
{
    /* The target precision is the number of bits in
     * double precision (53) + safety */
    const slong target_prec = 64;

    /* convert the input to arb_t
     * Since we are converting from binary (double precision)
     * to binary (arb_t), this conversion is exact
     */
    arb_t t_mp;
    arb_init(t_mp);
    arb_set_d(t_mp, t);
    assert(arb_is_exact(t_mp));

    arb_ptr F_mp = _arb_vec_init(m+1);

    slong working_prec = target_prec;
    int suff_acc = 0;

    while(!suff_acc)
    {
        working_prec += target_prec;
        suff_acc = 1;

        mirp_boys(F_mp, m, t_mp, working_prec);

        suff_acc = 1;
        for(int i = 0; i <= m; i++)
        {
            /* Do we have sufficient accuracy? We need at least
             * 53 bits + 11 bits safety OR the value is zero (has zero precision)
             * and the error bounds is exactly zero when converted to double precision */
            slong bits = arb_rel_accuracy_bits(F_mp + i);

            if(bits > 0 && bits < 64)
                suff_acc = 0;
            else if(bits <= 0)
            {
                double bound = mag_get_d(arb_radref(F_mp + i));
                if(bound > 0.0)
                    suff_acc = 0; 
            }
        }
    }

    /* convert back to double precision */
    for(int i = 0; i <= m; i++)
        F[i] = arf_get_d(arb_midref(F_mp + i), ARF_RND_NEAR);

    arb_clear(t_mp);
    _arb_vec_clear(F_mp, m+1);
}

