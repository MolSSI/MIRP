/*! \file
 *
 * \brief Miscellaneous helper functions for arblib data types
 */

#pragma once

#include <arb.h>

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Initialize an array of arb data types
 *
 * \param [inout] arr An array of <tt>arb_t</tt> to initialize
 * \param [in] length The number of elements in \p arr
 */
void mirp_init_arb_arr(arb_t * arr, int length);


/*! \brief Clear an array of arb data types
 *
 * \param [inout] arr An array of <tt>arb_t</tt> to initialize
 * \param [in] length The number of elements in \p arr
 */
void mirp_clear_arb_arr(arb_t * arr, int length);

#ifdef __cplusplus
}
#endif

