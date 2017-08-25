/*! \file
 *
 * \brief Calculation of the boys function using interval arithmetic
 */

#include "mirp/kernels/boys.h"
#include "mirp/math.h"
#include <arb_hypgeom.h>
#include <assert.h>

void mirp_boys(arb_ptr F, int m, const arb_t t, slong working_prec)
{
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

    /* et = exp(-x) */
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
        } while(!arb_contains(sum, test)); /* Is the old term contained completely
                                              within the new term */

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


int mirp_boys_target(arb_ptr F, int m, const arb_t t, slong target_prec)
{
    /* We run the boys function once, checking for the minimum
     * relative accuracy bits. If that is ok, we return.
     * If not, we enter a loop, increasing the working precision
     * until we reach our goal.
     *
     * We check if we enter an infinite loop, which can happen
     * if the precision of the inputs is not enough to
     * match the target prec. This is indicated by returning
     * nonzero
     */
    slong working_prec = target_prec + MIRP_BITS_INCREMENT;

    mirp_boys(F, m, t, working_prec);

    slong cur_bits = mirp_min_accuracy_bits(F, m+1);
    int suff_acc = (cur_bits >= target_prec);
    slong last_bits;

    while(!suff_acc)
    {
        last_bits = cur_bits;
        working_prec += MIRP_BITS_INCREMENT;

        mirp_boys(F, m, t, working_prec);

        cur_bits = mirp_min_accuracy_bits(F, m+1);

        if(cur_bits >= target_prec)
           suff_acc = 1;

        if(!suff_acc && cur_bits <= last_bits)
            break;
    }

    return !suff_acc;
}


void mirp_boys_target_str(arb_ptr F, int m, const char * t, slong target_prec)
{
    /* Procedure is similar to mirp_boys_target, however
     * this should always be guaranteed to succeed */
    arb_t t_mp;
    arb_init(t_mp);

    slong working_prec = target_prec;

    do
    {
        working_prec += MIRP_BITS_INCREMENT;
        arb_set_str(t_mp, t, working_prec);
        mirp_boys(F, m, t_mp, working_prec);

    } while(mirp_min_accuracy_bits(F, m+1) < target_prec);

    arb_clear(t_mp);
}


void mirp_boys_exact(double *F, int m, double t)
{
    /* The target precision is the number of bits in
     * double precision (53) + safety */
    const slong target_prec = 72;

    /* convert the input to arb_t
     * Since we are converting from binary (double precision)
     * to binary (arb_t), this conversion is exact
     */
    arb_t t_mp;
    arb_init(t_mp);
    arb_set_d(t_mp, t);
    assert(arb_is_exact(t_mp));

    arb_ptr F_mp = _arb_vec_init(m+1);

    mirp_boys_target(F_mp, m, t_mp, target_prec);

    /* convert back to double precision */
    for(int i = 0; i <= m; i++)
        F[i] = arf_get_d(arb_midref(F_mp + i), ARF_RND_NEAR);

    arb_clear(t_mp);
    _arb_vec_clear(F_mp, m+1);
}

