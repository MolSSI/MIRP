/*! \file
 *
 * \brief Some miscellaneous helper functions for mpfr data types
 */

#pragma once

#include <gmp.h>
#include <mpfr.h>

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Initialize an array of mpfr data types
 *
 * \param [inout] arr An array of <tt>mpfr_t</tt> to initialize
 * \param [in] length The number of elements in \p arr
 * \param [in] prec The precision to initialize the elements of \p arr with
 */
void mirp_init_mpfr_arr(mpfr_t * arr, int length, mpfr_prec_t prec);


/*! \brief Clear an array of mpfr data types
 *
 * \param [inout] arr An array of <tt>mpfr_t</tt> to initialize
 * \param [in] length The number of elements in \p arr
 */
void mirp_clear_mpfr_arr(mpfr_t * arr, int length);

#ifdef __cplusplus
}
#endif

