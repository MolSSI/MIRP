/*! \file
 *
 * \brief Sandbox for testing features
 */

#include "mirp_bin/cmdline.hpp"
#include <mirp/mirp.h>
#include <iostream>
#include <limits>

using namespace mirp;

/*! \brief Main function */
int main(void)
{
    double alpha1[3] = { 130.7093200, 23.8088610, 6.4436083 };
    double coeff1[3] = { 0.15432897, 0.53532814, 0.44463454 };
    double A[3] = { 0.0, 0.0, 0.0 };
    double integral[10000];

    mirp_eri_exact(integral,
                   3, A, 3, 1, alpha1, coeff1,
                   0, A, 3, 1, alpha1, coeff1,
                   0, A, 3, 1, alpha1, coeff1,
                   0, A, 3, 1, alpha1, coeff1);

    arb_ptr alpha_mp = _arb_vec_init(3);
    arb_ptr coeff_mp = _arb_vec_init(3);
    arb_ptr A_mp = _arb_vec_init(3);

    arb_ptr integral_mp = _arb_vec_init(10000);

    for(int i = 0; i < 3; i++)
    {
        arb_set_str(alpha_mp + i, std::to_string(alpha1[i]).c_str(), 48);
        arb_set_str(coeff_mp + i, std::to_string(coeff1[i]).c_str(), 48);
        arb_set_str(A_mp + i, std::to_string(A[i]).c_str(), 48);
    }

    int r = mirp_eri_target(integral_mp,
                            3, A_mp, 3, 1, alpha_mp, coeff_mp,
                            0, A_mp, 3, 1, alpha_mp, coeff_mp,
                            0, A_mp, 3, 1, alpha_mp, coeff_mp,
                            0, A_mp, 3, 1, alpha_mp, coeff_mp,
                            256);

    for(int i = 0; i < 10; i++)
    {
        printf("    Value %d: %32.18e\n", i, integral[i]);
        printf("Value ARB %d: [%d]  %s\n", i, r, arb_get_str(integral_mp + i, 100, 0));
    }

    return 0;
}
