#include "test_common.hpp"

#include <cmath>
#include <iostream>


namespace mirp {
namespace detail {

bool almost_equal(double a, double b, double tol)
{
    using std::fabs;
    using std::fmax;

    double diff = fabs(a-b);

    if(diff == 0)
        return true;

    double denominator = fmax(fabs(a), fabs(b));
    return (diff / denominator) < tol;
}


void print_results(unsigned long nfailed, unsigned long ntests)
{
    double percent_passed = 100.0 - (100.0 * (double)nfailed / (double)ntests);
    std::cout << nfailed << " / " << ntests << " failed ("
              << percent_passed << "% passed)\n";
}

} // close namespace detail
} // close namespace mirp
