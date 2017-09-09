/*! \file
 *
 * \brief mirp_run_test main function
 */

#include "mirp_bin/cmdline.hpp"
#include "mirp_bin/test_boys.hpp"
#include "mirp_bin/test_integral4.hpp"

#include <mirp/kernels/all.h>

#include <sstream>
#include <iostream>

using namespace mirp;

int main(int argc, char ** argv)
{
    arb_t a,b,c;
    arb_init(a);
    arb_init(b);
    arb_init(c);


    arb_set_str(a, "1.123456789012345678901234567890", 50);
    const char * a_s = arb_get_str(a, 5, ARB_STR_MORE);

    arb_set_str(b, a_s, 12);
    const char * b_s = arb_get_str(b, 5, ARB_STR_MORE);
    

    std::cout << "a: " << a_s << "\n";
    std::cout << "b: " << b_s << "\n";
    

    


    return 0;
}
