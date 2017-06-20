/*! \file
 *
 * \brief Some miscellaneous helper functions for mpfr data types
 */

#include "mirp/mpfr_help.h"
#include <assert.h>

void mirp_init_mpfr_arr(mpfr_t * arr, int length, mpfr_prec_t prec)
{
    assert(length >= 0);

    for(int i = 0; i < length; i++)
        mpfr_init2(arr[i], prec);
}


void mirp_clear_mpfr_arr(mpfr_t * arr, int length)
{
    assert(length >= 0);

    for(int i = 0; i < length; i++)
        mpfr_clear(arr[i]);
}

