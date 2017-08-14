/*! \file
 *
 * \brief Functions related to testing single ERI
 */

#pragma once

#include "mirp/typedefs.h"
#include <cmath>

namespace mirp {
namespace detail {


/*! \brief Compare two double-precision numbers for equality within a tolerance */
inline bool almost_equal(double a, double b, double tol)
{
    using std::fabs;
    using std::fmax;

    double diff = fabs(a-b);

    if(diff == 0)
        return true;

    double denominator = fmax(fabs(a), fabs(b));
    return (diff / denominator) < tol;
}


void print_results(unsigned long nfailed, unsigned long ntests);


void integral4_single_create_test(const std::string & input_filepath,
                                  const std::string & output_filepath,
                                  long ndigits, const std::string & header,
                                  cb_integral4_single_target_str cb);

long integral4_single_run_test(const std::string & filepath,
                               long target_prec,
                               cb_integral4_single_target_str cb);

long integral4_single_run_test_d(const std::string & filepath,
                                 cb_integral4_single_d cb);

long integral4_single_run_test_exact(const std::string & filepath,
                                     cb_integral4_single_exact cb,
                                     cb_integral4_single_target cb_mp);


} // close namespace detail
} // close namespace mirp

