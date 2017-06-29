/*! \file
 *
 * \brief Functions for manipulating shell structures
 */

#include "mirp/shell.h"

int mirp_iterate_gaussian(int lmn[3])
{
    const int am = lmn[0] + lmn[1] + lmn[2];
    if(lmn[2] >= am)
        return 0;

    if(lmn[2] < (am - lmn[0]))
    {
        lmn[1]--;
        lmn[2]++;
    }
    else
    {
        lmn[0]--;
        lmn[1] = am-lmn[0];
        lmn[2] = 0;
    }
    return 1;
}

#ifdef __cplusplus
}
#endif

