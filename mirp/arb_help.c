/*! \file
 *
 * \brief Some miscellaneous helper functions for arblib data types
 */

#include "mirp/arb_help.h"
#include <assert.h>

void mirp_init_arb_arr(arb_t * arr, int length)
{
    assert(length >= 0);

    for(int i = 0; i < length; i++)
        arb_init(arr[i]);
}


void mirp_clear_arb_arr(arb_t * arr, int length)
{
    assert(length >= 0);

    for(int i = 0; i < length; i++)
        arb_clear(arr[i]);
}

