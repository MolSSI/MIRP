#include "test_common.hpp"

#include <iostream>


namespace mirp {
namespace detail {

void print_results(unsigned long nfailed, unsigned long ntests)
{
    double percent_passed = 100.0 - (100.0 * (double)nfailed / (double)ntests);
    std::cout << nfailed << " / " << ntests << " failed ("
              << percent_passed << "% passed)\n";
}

} // close namespace detail
} // close namespace mirp
