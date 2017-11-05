/*! \file
 *
 * \brief Aggregate header for callback helpers
 */

#pragma once

namespace mirp {


/*! \brief Helps with calling callbacks with the given number of centers
 *
 * This structure is used to appropriately unpack arrays of data (of length \p N)
 * when calling C-style callbacks.
 *
 * \tparam Number of centers
 */
template<int N> struct callback_helper;

} // close namespace mirp

#include "mirp_bin/callback_helper4.hpp"
